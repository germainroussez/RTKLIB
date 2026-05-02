/*------------------------------------------------------------------------------
* test_has_mt1.c : standalone validation of the HAS MT1 message parser
*
*          Copyright (C) 2026 by G.ROUSSEZ, All rights reserved.
*
* description :
*     Three checks are run:
*
*       (1) A minimal MT1 carrying only a Mask block (one Galileo GNSS, two
*           satellites, two signals, no cell mask) is hand-packed bit by bit
*           and parsed back. Every header and Mask field must round-trip.
*
*       (2) A complete MT1 with Mask, Orbit Corrections, Clock Full-Set
*           Corrections, Code Biases and Phase Biases is hand-packed for two
*           Galileo satellites and two signals, with the cell mask enabled
*           and a couple of cells masked off, and parsed back. Every value
*           must match the input bit for bit.
*
*       (3) Smoke test: the 795-octet message produced by the Phase 3
*           pipeline (ICD Annex C reference example) is fed to
*           has_parse_mt1(). The parser must accept the message, recognise
*           the MT1 header, walk the blocks without truncation, and yield
*           plausible field values.
*
* references :
*     [1] European GNSS Service Centre, Galileo High Accuracy Service
*         Signal-in-Space Interface Control Document, Issue 1.0, May 2022.
*
* usage :
*     cc -O2 -I../../src ../../src/has.c ../../src/rtkcmn.c \
*        test_has_mt1.c -o test_has_mt1
*     ./test_has_mt1
*
* version : $Revision:$ $Date:$
* history : 2026/05/02 1.0 new
*-----------------------------------------------------------------------------*/
#include "rtklib.h"
#include <stdio.h>
#include <string.h>

/* small helper to pack a bit-stream incrementally ---------------------------*/
typedef struct {
    uint8_t buf[4096];
    int off;            /* current bit offset */
} packer_t;

static void pack_u(packer_t *p, int nbits, uint32_t v)
{
    setbitu(p->buf, p->off, nbits, v);
    p->off += nbits;
}
static void pack_s(packer_t *p, int nbits, int32_t v)
{
    /* mask off the high bits, getbitu/sign_extend will recover the value    */
    uint32_t u = (uint32_t)v & ((nbits == 32) ? 0xFFFFFFFFu :
                                ((1u << nbits) - 1));
    setbitu(p->buf, p->off, nbits, u);
    p->off += nbits;
}
/* pack a 40-bit satellite mask with the satellites listed in svids[]        */
/* (1-based per ICD Table 19); MSB of the field corresponds to SVID=1.       */
static void pack_sat_mask(packer_t *p, const int *svids, int n)
{
    uint8_t m[5] = {0};
    int i;
    for (i = 0; i < n; i++) {
        int svid = svids[i];
        if (svid >= 1 && svid <= 40) setbitu(m, svid - 1, 1, 1);
    }
    for (i = 0; i < 5; i++) pack_u(p, 8, m[i]);
}
/* pack a 16-bit signal mask from a list of signal indexes (0-based, ICD 5.2.1.3)*/
static void pack_sig_mask(packer_t *p, const int *sigs, int n)
{
    uint16_t m = 0;
    int i;
    for (i = 0; i < n; i++) {
        int s = sigs[i];
        if (s >= 0 && s <= 15) m |= (uint16_t)(1u << (15 - s));
    }
    pack_u(p, 16, m);
}
/* Test 1: minimal MT1 with just the Mask block ------------------------------*/
static int test_minimal_mask(void)
{
    packer_t p = {{0}, 0};
    gal_has_msg_t out;
    int svids[] = {3, 18};
    int sigs[]  = {0, 8};               /* E1-B + E5b-I+E5b-Q */
    int n;

    /* MT1 header: TOH=42, mask=1, others=0, mask_id=7, iod=11 -------------- */
    pack_u(&p, 12, 42);
    pack_u(&p, 1, 1);                   /* mask flag */
    pack_u(&p, 1, 0);
    pack_u(&p, 1, 0);
    pack_u(&p, 1, 0);
    pack_u(&p, 1, 0);
    pack_u(&p, 1, 0);
    pack_u(&p, 4, 0);                   /* reserved */
    pack_u(&p, 5, 7);                   /* mask_id */
    pack_u(&p, 5, 11);                  /* iod_set_id */

    /* Mask block: Nsys=1, Galileo, 2 sats, 2 sigs, no cell mask, NM=0 ------ */
    pack_u(&p, 4, 1);                   /* nsys */
    pack_u(&p, 4, HAS_GNSS_GAL);
    pack_sat_mask(&p, svids, 2);
    pack_sig_mask(&p, sigs, 2);
    pack_u(&p, 1, 0);                   /* cmaf = 0 */
    pack_u(&p, 3, 0);                   /* nm = 0 */
    pack_u(&p, 6, 0);                   /* mask block reserved trailer */

    n = (p.off + 7) / 8;
    if (has_parse_mt1(p.buf, n, &out) != 0) {
        fprintf(stderr, "test1: has_parse_mt1 returned non-zero\n");
        return 1;
    }
    if (out.hdr.toh != 42 || out.hdr.mask_id != 7 || out.hdr.iod_set_id != 11
            || !out.hdr.mask_flag || out.hdr.orbit_flag) {
        fprintf(stderr, "test1: header round-trip failed\n");
        return 1;
    }
    if (out.nsys != 1 || out.mask[0].gnss_id != HAS_GNSS_GAL
            || out.mask[0].nsat != 2 || out.mask[0].nsig != 2
            || out.mask[0].sat_idx[0] != 3 || out.mask[0].sat_idx[1] != 18
            || out.mask[0].sig_idx[0] != 0 || out.mask[0].sig_idx[1] != 8
            || out.mask[0].cmaf != 0 || out.mask[0].nm != 0) {
        fprintf(stderr,
            "test1: mask round-trip failed (gnss=%d nsat=%d nsig=%d "
            "sat=[%d,%d] sig=[%d,%d] cmaf=%d nm=%d)\n",
            out.mask[0].gnss_id, out.mask[0].nsat, out.mask[0].nsig,
            out.mask[0].sat_idx[0], out.mask[0].sat_idx[1],
            out.mask[0].sig_idx[0], out.mask[0].sig_idx[1],
            out.mask[0].cmaf, out.mask[0].nm);
        return 1;
    }
    printf("test1: minimal MT1 mask round-trip OK\n");
    return 0;
}
/* Test 2: full MT1 with mask + orbit + clock + code bias + phase bias ------*/
static int test_complete_mt1(void)
{
    packer_t p = {{0}, 0};
    gal_has_msg_t out;
    int svids[] = {5, 27};
    int sigs[]  = {0, 5};                  /* E1-B + E5a-I+E5a-Q */
    int n;
    /* synthetic values to round-trip */
    int orbit_iod[] = {123, 456};
    int orbit_dr[]  = { 100, -200};        /* raw, scale 0.0025 m */
    int orbit_dit[] = { -50,  500};        /* raw, scale 0.008 m  */
    int orbit_dct[] = {  10,  -10};        /* raw, scale 0.008 m  */
    int dcm = 2;                           /* multiplier 1..4, raw 0..3 */
    int dcc[]      = { 300,  -75};         /* raw, scale 0.0025 m */
    int cb_vals[][2] = { { 50, -100}, {  25,    0} };  /* sat x sig, raw 0.02m */
    int pb_vals[][2] = { {-30,  120}, { 200, -200} };  /* sat x sig, raw 0.01c */
    int pdi_vals[][2]= { {  1,    2}, {   3,    0} };
    int cell_mask_bits[][2] = { {1, 1}, {1, 0} };       /* drop one cell */

    /* MT1 header: all flags except clk_sub */
    pack_u(&p, 12, 1234);
    pack_u(&p, 1, 1);   /* mask */
    pack_u(&p, 1, 1);   /* orbit */
    pack_u(&p, 1, 1);   /* clk full */
    pack_u(&p, 1, 0);   /* clk sub */
    pack_u(&p, 1, 1);   /* cb */
    pack_u(&p, 1, 1);   /* pb */
    pack_u(&p, 4, 0);
    pack_u(&p, 5, 17);  /* mask_id */
    pack_u(&p, 5, 25);  /* iod_set_id */

    /* Mask block: Galileo + cmaf=1 ----------------------------------------- */
    pack_u(&p, 4, 1);
    pack_u(&p, 4, HAS_GNSS_GAL);
    pack_sat_mask(&p, svids, 2);
    pack_sig_mask(&p, sigs, 2);
    pack_u(&p, 1, 1);                      /* cmaf */
    {
        int row, col;
        for (row = 0; row < 2; row++) {
            for (col = 0; col < 2; col++) {
                pack_u(&p, 1, cell_mask_bits[row][col]);
            }
        }
    }
    pack_u(&p, 3, 0);                      /* nm */
    pack_u(&p, 6, 0);                      /* reserved */

    /* Orbit block: VI=5 (60s) + 2 SV ---------------------------------------*/
    pack_u(&p, 4, 5);
    {
        int i;
        for (i = 0; i < 2; i++) {
            pack_u(&p, 10, orbit_iod[i]);  /* GAL IOD = 10 bits */
            pack_s(&p, 13, orbit_dr[i]);
            pack_s(&p, 12, orbit_dit[i]);
            pack_s(&p, 12, orbit_dct[i]);
        }
    }
    /* Clock Full-Set block: VI=5, DCM=2, 2 DCC -----------------------------*/
    pack_u(&p, 4, 5);
    pack_u(&p, 2, dcm - 1);                /* DCM raw 0..3 */
    {
        int i;
        for (i = 0; i < 2; i++) pack_s(&p, 13, dcc[i]);
    }
    /* Code Bias block: VI=4, only cells where cell_mask_bits == 1 ---------*/
    pack_u(&p, 4, 4);
    {
        int row, col;
        for (row = 0; row < 2; row++) {
            for (col = 0; col < 2; col++) {
                if (cell_mask_bits[row][col]) {
                    pack_s(&p, 11, cb_vals[row][col]);
                }
            }
        }
    }
    /* Phase Bias block: VI=4, same cell mask --------------------------------*/
    pack_u(&p, 4, 4);
    {
        int row, col;
        for (row = 0; row < 2; row++) {
            for (col = 0; col < 2; col++) {
                if (cell_mask_bits[row][col]) {
                    pack_s(&p, 11, pb_vals[row][col]);
                    pack_u(&p, 2, pdi_vals[row][col]);
                }
            }
        }
    }
    n = (p.off + 7) / 8;

    if (has_parse_mt1(p.buf, n, &out) != 0) {
        fprintf(stderr, "test2: has_parse_mt1 returned non-zero\n");
        return 1;
    }
    if (out.hdr.toh != 1234 || out.hdr.mask_id != 17 ||
        out.hdr.iod_set_id != 25) {
        fprintf(stderr, "test2: header values wrong\n");
        return 1;
    }
    if (out.orbit_vi != 5 || out.clk_full_vi != 5 ||
        out.cb_vi != 4 || out.pb_vi != 4) {
        fprintf(stderr, "test2: VI values wrong\n");
        return 1;
    }
    if (out.dcm[0] != dcm) {
        fprintf(stderr, "test2: DCM=%d expected %d\n", out.dcm[0], dcm);
        return 1;
    }
    {
        int i;
        for (i = 0; i < 2; i++) {
            if (out.orbit[0][i].iod != orbit_iod[i] ||
                out.orbit[0][i].dr  != orbit_dr[i]  ||
                out.orbit[0][i].dit != orbit_dit[i] ||
                out.orbit[0][i].dct != orbit_dct[i] ||
                out.dcc[0][i] != dcc[i]) {
                fprintf(stderr, "test2: orbit/clock SV%d round-trip failed\n",
                        i);
                return 1;
            }
        }
    }
    {
        int row, col;
        for (row = 0; row < 2; row++) {
            for (col = 0; col < 2; col++) {
                int valid = cell_mask_bits[row][col];
                const gal_has_bias_t *b = &out.bias[0][row][col];
                if (valid && (b->cb != cb_vals[row][col] ||
                              !b->cb_valid)) {
                    fprintf(stderr, "test2: CB row=%d col=%d wrong\n",
                            row, col);
                    return 1;
                }
                if (valid && (b->pb != pb_vals[row][col] ||
                              b->pdi != pdi_vals[row][col] ||
                              !b->pb_valid)) {
                    fprintf(stderr, "test2: PB row=%d col=%d wrong\n",
                            row, col);
                    return 1;
                }
                if (!valid && (b->cb_valid || b->pb_valid)) {
                    fprintf(stderr,
                        "test2: cell-masked-off row=%d col=%d still flagged\n",
                        row, col);
                    return 1;
                }
            }
        }
    }
    printf("test2: complete MT1 round-trip OK (mask + orbit + clock + "
           "biases, %d bits)\n", p.off);
    return 0;
}
/* The 795-octet message recovered by the Phase 3 pipeline from ICD Annex C. */
/* This is the byte-for-byte output of has_rs_decode() applied to the 15     */
/* encoded pages of Annex C. Used here for a syntactic smoke test only.      */
static const uint8_t annexc_msg[795] = {
  0, 12,192, 11, 32,255,223,255,255,  0,129,  0,247,255,255,125,245, 95,253,254,
 11,238,232,167,154, 65, 36, 16,  0,166,  0, 10,  1,160, 18,128, 64,  2,  0, 32,
  1, 19,251,192, 65,254,187,240,  0,128,  8,  0, 66,
255,104, 34,254,162, 24,  7,193,147,247, 89,128, 53,253,127,106, 47,  0,  8,  0,
128,  1,111,249,  2,135,231,150,127,112, 37,128, 88,127,238, 33,122, 16,201,223,
204, 14,127,101, 29,245,119,217,129, 96, 63,254, 65,
 71,249,  3,255,157,247,128, 92, 21,255,159,220,255,128,  8,  0, 64,  4,  0, 10,
  0, 36,  7,255,157,124,  7,223,127,254, 43, 95,220,238, 48, 85, 25,  1, 31,215,
253, 36, 71,159,  0, 80, 14,142,126,220, 49, 64, 28,
 67,253,176, 35,  4,  0,127,229,  3, 15,241,172, 64,  2,  0, 32,  0,  2,  0, 16,
  1,  0,  7,127,236,  6,224,  1, 65,254,176, 42,252,178,196,  0, 32,  2,  0,  4,
 63,245,246,192, 34,  9,127,124, 14, 63, 68, 18,255,
 79,225,255,136, 37,254,143,252,255,  0, 72,  8, 31,227,253,160,151,244,192, 75,
243,129, 47,229,255, 39,240,  2, 95,198,255, 95,244,  4,128,237,250, 96, 28,  8,
255,232,  2, 63,204, 15,  0,176, 11,128,168, 37,253,
240, 15,255,112, 75,247, 31,255,253,192,151,251, 64, 12,  0,129, 47,231,129,167,
248,  2, 95,230,  2, 32, 50,  4,128, 16,  1,160, 22,  7,255,208,  6, 64, 64, 18,
254,192, 14,  0,  8, 37,252,127,229,  0,192, 75,255,
 64, 86,  5,192,136,  4,  0, 68,  3,  1, 47,226,127,239,251,240,187, 35,220,148,
 69,142,240, 66, 10,254, 31,166, 21, 68,171,218,119,193, 48, 68, 67, 32,161, 16,
 67,  3,211,247,111,101,251,190,231,204,245,254,107,
221,248,191,207,244,121,183,165,241,220, 59,243,252,225, 36, 59, 68,233, 13, 23,
132,172, 53, 11, 47, 41,242,189, 96,123, 26, 30,123,178,  7, 81,146,  1,  0, 56,
  7,  6,159,143,235,124,240, 12, 13, 66,216, 91,  6,
 31, 51,210,250,127,160, 15,195, 80,106,  2,  1, 92, 75,  9, 64,155,240,124,191,
149,  4,  0,100, 21,130,160, 79,200,244, 14,136,210,221,159,115,239,189,196,  0,
128, 64,  4,  7,193,152, 88,138,208,233,244, 61,103,
174,249,  0,156, 34,  4, 32,205,239,188,159,144,249, 32,240, 51,134, 96, 64, 26,
 69,160,180, 17,160,132, 28,131,128,194,  6,193,136, 45,  1, 33, 36, 62,135,208,
 43,242,125, 31,162,252, 97,132, 81,138, 80,220,176,
  0,128,  4,  0, 32,  1,  0,  8,  0, 64,  2,  0, 16,  0,128,  4,  0, 32,  1,  0,
  8,  0, 64,  2,  0, 16,  0,128,  4,  0, 32,  1,  0,  8,  0, 64,  2,  0, 16,  0,
128,  4,  0, 32,  1,  0,  8,  0, 64,  2,  0, 16,  0,
128,  4,  0, 32,  1,  0,  8,  0, 64,  2,  0, 16,  0,128,  4,  0, 32,  1,  0,  8,
  0, 64,  2,  0, 16,  0,128,  4,  0, 32,  1,  0,  8,  0, 64,  2,  0, 16,  0,128,
  4,  0, 32,  1,  0,  8,  0, 64,  2,  0, 16,  0,128,
  4,  0, 32,  1,  0,  8,  0, 64,  2,  0, 16,  0,128,  4,  0, 32,  1,  0,  8,  0,
 64,  2,  0, 16,  0,128,  4,  0, 32,  1,  0,  8,  0, 64,  2,  0, 16,  0,128,  4,
  0, 32,  1,  0,  8,  0, 64,  2,  0, 16,  0,128,  4,
  0, 32,  1,  0,  8,  0, 64,  2,  0, 16,  0,128,  4,  0, 32,  1,  0,  8,  0, 64,
  2,  0, 16,  0,128,  4,  0, 32,  1,  0,  8,  0, 64,  2,  0, 16,  0,128,  4,  0,
 32,  1,  0,  8,  0, 64,  2,  0, 16,  0,128,  4,  0,
 32,  1,  0,  8,  0, 64,  2,  0, 16,  0,128,  4,  0, 32,  1,  0,  8,  0, 42,170,
170,170,170,170,170,170,170,170,170,170,170,170,170,170,170,170,170,170,170,170,
170,170,170,170,170,170,170,170,170,170,170,170,170
};

static int test_annexc_smoke(void)
{
    gal_has_msg_t out;
    int rc;

    rc = has_parse_mt1(annexc_msg, (int)sizeof(annexc_msg), &out);
    if (rc != 0) {
        fprintf(stderr, "test3: parse_mt1 returned %d\n", rc);
        return 1;
    }
    /* The Annex C header bytes 00 0C C0 0B decode to:                       */
    /*   TOH=0, mask=1 orbit=1 clk_full=0 clk_sub=0 cb=1 pb=1, mask_id=0,    */
    /*   iod_set_id=11.                                                       */
    if (out.hdr.toh != 0 || out.hdr.mask_id != 0 || out.hdr.iod_set_id != 11
            || !out.hdr.mask_flag || !out.hdr.orbit_flag
            || out.hdr.clk_full_flag || out.hdr.clk_sub_flag
            || !out.hdr.cb_flag || !out.hdr.pb_flag) {
        fprintf(stderr,
            "test3: unexpected MT1 header (toh=%u mask=%u orbit=%u "
            "clk_full=%u clk_sub=%u cb=%u pb=%u mask_id=%u iod=%u)\n",
            out.hdr.toh, out.hdr.mask_flag, out.hdr.orbit_flag,
            out.hdr.clk_full_flag, out.hdr.clk_sub_flag,
            out.hdr.cb_flag, out.hdr.pb_flag,
            out.hdr.mask_id, out.hdr.iod_set_id);
        return 1;
    }
    if (out.nsys < 1 || out.nsys > HAS_NSYS_MAX) {
        fprintf(stderr, "test3: implausible nsys=%d\n", out.nsys);
        return 1;
    }
    printf("test3: Annex C 795-octet message parsed OK "
           "(nsys=%d mask0 nsat=%d nsig=%d)\n",
           out.nsys, out.mask[0].nsat, out.mask[0].nsig);
    return 0;
}
int main(void)
{
    int rc = 0;
    rc |= test_minimal_mask();
    rc |= test_complete_mt1();
    rc |= test_annexc_smoke();
    if (rc) {
        fprintf(stderr, "test_has_mt1: FAIL\n");
        return 1;
    }
    printf("test_has_mt1: ALL TESTS PASS\n");
    return 0;
}
