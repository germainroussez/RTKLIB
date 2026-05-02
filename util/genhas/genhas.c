/*------------------------------------------------------------------------------
* genhas.c : generate Galileo HAS Reed-Solomon constants
*
*          Copyright (C) 2026 by G.ROUSSEZ, All rights reserved.
*
* description :
*     This standalone tool generates the constant tables used by the
*     Reed-Solomon (255,32,224) erasure decoder for the Galileo HAS High
*     Parity Vertical Reed-Solomon (HPVRS) outer-layer coding scheme. The
*     output is C source intended to be pasted into src/has.c.
*
*     Three tables are produced:
*       gf_exp[512]  : alpha^i for i=0..510 over GF(256) defined by
*                      p(alpha) = alpha^8 + alpha^4 + alpha^3 + alpha^2 + 1
*                      (octet 0x1D, leading bit implied at bit 8). Doubled
*                      to 512 entries to allow gf_mul without modulo 255.
*       gf_log[256]  : log_alpha(x) for x=1..255, with gf_log[0] unused.
*       has_g[255*32]: row-major systematic generator matrix G of the
*                      RS(255,32) code, in the indexing convention of
*                      sections 6.2.3 and 6.3 of the HAS SIS ICD; G is
*                      255 rows by 32 columns of GF(256) octets, with
*                      the upper 32 rows forming an identity submatrix
*                      and the lower 223 rows forming the parity submatrix.
*
* references :
*     [1] European GNSS Service Centre, Galileo High Accuracy Service Signal-
*         in-Space Interface Control Document, Issue 1.0, May 2022.
*
* usage :
*     cc -O2 -o genhas genhas.c
*     ./genhas > has_rs_tables.c
*
* version : $Revision:$ $Date:$
* history : 2026/05/02 1.0 new
*-----------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdint.h>
#include <string.h>

/* primitive polynomial p(alpha) = alpha^8+alpha^4+alpha^3+alpha^2+1, ICD Eq.7 */
#define HAS_RS_PRIM     0x11D

/* RS code parameters, ICD section 6.2 ---------------------------------------*/
#define HAS_RS_N        255
#define HAS_RS_K         32

static uint8_t gf_exp[512];
static uint8_t gf_log[256];
static uint8_t has_g[HAS_RS_N * HAS_RS_K];

/* GF(256) multiplication using log/exp tables -------------------------------*/
static uint8_t gf_mul(uint8_t a, uint8_t b)
{
    if (a == 0 || b == 0) return 0;
    return gf_exp[gf_log[a] + gf_log[b]];
}
/* build gf_exp[] / gf_log[] under the HAS primitive polynomial --------------*/
static void build_gf_tables(void)
{
    int i;
    uint16_t x = 1;

    for (i = 0; i < 255; i++) {
        gf_exp[i] = (uint8_t)x;
        gf_log[(uint8_t)x] = (uint8_t)i;
        x <<= 1;
        if (x & 0x100) x ^= HAS_RS_PRIM;
    }
    /* duplicate exp table so callers can sum two logs without modulo */
    for (i = 255; i < 512; i++) gf_exp[i] = gf_exp[i - 255];
}
/* systematic RS encoding of a length-K information vector -------------------*/
/* Encodes c = [c_31, c_30, ..., c_0]^T into Gamma = [Gamma_254, ..., Gamma_0] */
/* per ICD Eq. 11/12: Gamma = (c, gamma) with gamma the parity from r(x) =    */
/* c(x)*x^{n-k} mod g(x). The K information rows of Gamma are exactly c, in   */
/* the index ordering that yields G as a (I_K stacked over P_{N-K}) systematic*/
/* generator. We compute g(x) on the fly as g(x) = prod_{i=1..N-K}(x-alpha^i) */
static void rs_encode_systematic(const uint8_t *info_K, uint8_t *codeword_N)
{
    /* generator polynomial g(x), degree N-K = 223, coeffs in ascending x */
    static uint8_t gpoly[HAS_RS_N - HAS_RS_K + 1];
    static int gpoly_built = 0;
    int i, j;

    if (!gpoly_built) {
        memset(gpoly, 0, sizeof(gpoly));
        gpoly[0] = 1;
        for (i = 1; i <= HAS_RS_N - HAS_RS_K; i++) {
            /* multiply current gpoly by (x - alpha^i): in GF(256) -alpha^i = alpha^i */
            uint8_t alpha_i = gf_exp[i];
            for (j = i; j > 0; j--) {
                gpoly[j] = gpoly[j - 1] ^ gf_mul(gpoly[j], alpha_i);
            }
            gpoly[0] = gf_mul(gpoly[0], alpha_i);
        }
        gpoly_built = 1;
    }
    /* parity = (c(x) * x^{N-K}) mod g(x), implemented as long division */
    static uint8_t parity[HAS_RS_N - HAS_RS_K];
    memset(parity, 0, sizeof(parity));

    /* c(x) = sum c_j * x^j with j = 0..K-1 ; info_K[i] = c_{K-1-i} (i.e.   */
    /* info_K[0] = c_{K-1} = highest power coefficient first per Eq.10).    */
    /* We feed the coefficients high-to-low into the LFSR-style division.   */
    for (i = 0; i < HAS_RS_K; i++) {
        uint8_t feedback = info_K[i] ^ parity[HAS_RS_N - HAS_RS_K - 1];
        for (j = HAS_RS_N - HAS_RS_K - 1; j > 0; j--) {
            parity[j] = parity[j - 1] ^ gf_mul(feedback, gpoly[j]);
        }
        parity[0] = gf_mul(feedback, gpoly[0]);
    }
    /* assemble codeword in ICD index ordering Gamma_254 down to Gamma_0:    */
    /* first K entries = c (information part), next N-K entries = parity     */
    /* in descending power order (gamma_{N-K-1}, ..., gamma_0).              */
    for (i = 0; i < HAS_RS_K; i++) codeword_N[i] = info_K[i];
    for (i = 0; i < HAS_RS_N - HAS_RS_K; i++) {
        codeword_N[HAS_RS_K + i] = parity[HAS_RS_N - HAS_RS_K - 1 - i];
    }
}
/* build the full 255 x 32 systematic generator matrix G ---------------------*/
/* G is obtained by encoding each of the K canonical unit vectors e_j and    */
/* taking the resulting N-length codewords as the columns of G; equivalently */
/* the upper 32 x 32 submatrix is the identity, and the lower 223 x 32 block */
/* P is the parity contribution. Stored row-major: has_g[r*32 + c].          */
static void build_generator_matrix(void)
{
    int j, r;
    uint8_t info[HAS_RS_K];
    uint8_t code[HAS_RS_N];

    for (j = 0; j < HAS_RS_K; j++) {
        memset(info, 0, sizeof(info));
        info[j] = 1;            /* unit vector with 1 at row j */
        rs_encode_systematic(info, code);
        for (r = 0; r < HAS_RS_N; r++) {
            has_g[r * HAS_RS_K + j] = code[r];
        }
    }
}
/* print a uint8_t table in groups of 16 per line ----------------------------*/
static void dump_u8_table(const char *name, const uint8_t *tbl, int n,
                          const char *comment)
{
    int i;
    printf("/* %s */\n", comment);
    printf("static const uint8_t %s[%d] = {\n", name, n);
    for (i = 0; i < n; i++) {
        if ((i % 16) == 0) printf("    ");
        printf("%3u%s", tbl[i], i == n - 1 ? "" : ",");
        if ((i % 16) == 15 || i == n - 1) printf("\n");
        else printf(" ");
    }
    printf("};\n");
}
/* dump the 255 x 32 generator matrix as a flat row-major table --------------*/
static void dump_g_matrix(void)
{
    int r, c;
    printf("/* HAS RS(255,32) systematic generator matrix G, row-major.\n");
    printf(" * Row r in 0..254 corresponds to encoded page PID = r+1.        \n");
    printf(" * Column c in 0..31 corresponds to information symbol c_{K-1-c}.\n");
    printf(" * Rows 0..31 form the identity submatrix per ICD section 6.2.3. */\n");
    printf("static const uint8_t has_rs_g[%d][%d] = {\n", HAS_RS_N, HAS_RS_K);
    for (r = 0; r < HAS_RS_N; r++) {
        printf("    {");
        for (c = 0; c < HAS_RS_K; c++) {
            printf("%3u%s", has_g[r * HAS_RS_K + c],
                   c == HAS_RS_K - 1 ? "" : ",");
            if (c != HAS_RS_K - 1) printf(" ");
        }
        printf("}%s\n", r == HAS_RS_N - 1 ? "" : ",");
    }
    printf("};\n");
}
/* sanity check against ICD Annex C decoding example -------------------------*/
/* The ICD provides the 15 x 15 sub-matrix D obtained by selecting K=MS=15   */
/* received PIDs and the first 15 columns of G. We rebuild D from has_g and  */
/* compare element-wise to the values printed in the ICD.                   */
static int annex_c_self_check(void)
{
    static const int pids[15] = {55,56,57,58,59,174,175,176,187,188,
                                 239,240,241,252,253};
    static const uint8_t expected[15][15] = {
        { 31, 50,155,253,213,220, 84,174,239, 85, 87,105,214, 81,160},
        {113, 18, 35,135,205, 43,156, 23,127,169,162,160, 15, 49,202},
        {204,239,127,208, 89,187, 30,192, 37,152,221,214,211, 49, 93},
        { 72,  7, 24, 67,  1,245,154,234, 84,179, 37, 96,222, 33, 64},
        {253,151,182,118,101,136,118,241,195, 26,152, 14,225, 28,193},
        {114,171,242,238, 47,124, 59,125, 65, 23, 39,150,161,226,  5},
        { 33, 32,  3,  8, 36,151,121, 17,218, 26, 98, 82, 65,146,162},
        { 37,190,149, 41, 64, 68,119, 19,153, 51,235,147,203,136,225},
        { 58,217, 47, 14,  1, 13,117,  8,167, 10,105,226, 96,158,229},
        {169,119,204,119, 80, 22, 46, 55,120, 70, 39, 68,156,140,150},
        {145, 19,150, 65,190, 97,199,178, 76,115,138,198,136, 18,180},
        {235,120, 75, 39,150,196, 72,209,145, 27,180, 77, 11,  2,154},
        {143,165, 24,101,222,187,133, 80,114, 98,164, 11, 16,227, 43},
        { 15,105,201,161,101,197,235,191,127, 28,238,232,231,198,234},
        { 84,157,205,255,217,251,101,194,230,208, 26,232, 23,201, 46}
    };
    int i, j, mismatches = 0;

    for (i = 0; i < 15; i++) {
        for (j = 0; j < 15; j++) {
            uint8_t got = has_g[(pids[i] - 1) * HAS_RS_K + j];
            if (got != expected[i][j]) {
                fprintf(stderr,
                    "MISMATCH at PID %d, col %d: expected %u, got %u\n",
                    pids[i], j, expected[i][j], got);
                mismatches++;
            }
        }
    }
    return mismatches;
}
int main(int argc, char **argv)
{
    int check_only = (argc > 1 && strcmp(argv[1], "-c") == 0);
    int mismatches;

    build_gf_tables();
    build_generator_matrix();

    mismatches = annex_c_self_check();
    if (mismatches != 0) {
        fprintf(stderr,
            "genhas: ICD Annex C cross-check FAILED (%d mismatches)\n",
            mismatches);
        return 2;
    }
    fprintf(stderr,
        "genhas: ICD Annex C cross-check OK (15x15 D-matrix matches)\n");

    if (check_only) return 0;

    printf("/* Galileo HAS RS(255,32,224) precomputed constants. */\n");
    printf("/* Generated by util/genhas/genhas.c - DO NOT EDIT BY HAND. */\n");
    printf("/* Reference: Galileo HAS SIS ICD v1.0, May 2022, section 6. */\n\n");
    dump_u8_table("has_rs_gf_exp", gf_exp, 512,
                  "alpha^i over GF(256), i=0..510, doubled to allow no-mod mul");
    printf("\n");
    dump_u8_table("has_rs_gf_log", gf_log, 256,
                  "log_alpha(i) for i=1..255, gf_log[0] unused");
    printf("\n");
    dump_g_matrix();
    return 0;
}
