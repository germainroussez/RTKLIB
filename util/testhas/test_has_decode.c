/*------------------------------------------------------------------------------
* test_has_decode.c : standalone validation of the HAS page parser and
*                     multi-page reassembly state machine in src/has.c
*
*          Copyright (C) 2026 by G.ROUSSEZ, All rights reserved.
*
* description :
*     Three checks are run:
*
*       (1) The C/NAV page printed in HAS SIS ICD Annex C is shifted to the
*           transport-agnostic 56-octet HAS Page form expected by
*           has_input_page() and fed to it. The parser must identify it as a
*           non-dummy page with the expected header values (HASS=0, MT=1,
*           MID=15, MS=15, PID=55) and store its encoded payload at index 0
*           of channel 0.
*
*       (2) A 56-octet HAS Dummy Page (first three bytes = 0xAF, 0x3B, 0xC3)
*           is fed to has_input_page() and must be silently discarded.
*
*       (3) A complete 15-page reassembly is exercised end-to-end. The 15
*           encoded payloads from ICD Annex C are wrapped into synthetic
*           HAS Pages (header MT=1 MID=7 MS=15 with the Annex C PID list,
*           encoded payload appended after the 24-bit header), pushed in
*           shuffled order, and the resulting message must match the
*           15-page reference output of Annex C byte for byte.
*
* references :
*     [1] European GNSS Service Centre, Galileo High Accuracy Service
*         Signal-in-Space Interface Control Document, Issue 1.0, May 2022,
*         Annex C (Reed-Solomon Decoding Example).
*
* usage :
*     cc -O2 -I../../src ../../src/has.c ../../src/rtkcmn.c \
*        test_has_decode.c -o test_has_decode
*     ./test_has_decode
*
* version : $Revision:$ $Date:$
* history : 2026/05/02 1.0 new
*-----------------------------------------------------------------------------*/
#include "rtklib.h"
#include <stdio.h>
#include <string.h>

/* ICD Annex C reference C/NAV page (64 octets, raw u-blox SFRBX form).      */
/* The HAS Page is at bit offset 14 within this buffer; we shift it down to  */
/* 56 byte-aligned octets in cnav_to_has_page() before feeding has_input_page. */
static const uint8_t annexc_cnav[64] = {
    0xff,0xfc,0x17,0xb8, 0xde,0x11,0xef,0x1d,
    0x27,0xad,0xf5,0xc5, 0xd0,0x91,0x1e,0x23,
    0xed,0x15,0x1a,0x46, 0x30,0x00,0x9c,0xab,
    0xaf,0x05,0x52,0x4b, 0x31,0xba,0xd5,0x69,
    0x62,0x03,0x89,0x86, 0xeb,0x8c,0x5c,0x68,
    0x8f,0x74,0x2f,0x95, 0x8b,0xf2,0x35,0xbf,
    0x62,0x39,0x88,0xa7, 0x0a,0x79,0xf6,0x32,
    0x67,0x7d,0x0c,0x46, 0x90,0x00,0x00,0x00
};

/* The 15 received PIDs and encoded payloads from ICD Annex C, used in test 3*/
static const uint8_t pids_anc[15] = {
    55, 56, 57, 58, 59, 174, 175, 176, 187, 188, 239, 240, 241, 252, 253
};
static const uint8_t enc_anc[15][53] = {
{132,123,199, 73,235,125,113,116, 36, 71,136,251, 69, 70,145,140,  0, 39, 42,235,
 193, 84,146,204,110,181, 90, 88,128,226, 97,186,227, 23, 26, 35,221, 11,229, 98,
 252,141,111,216,142, 98, 41,194,158,125,140,153,223},
{ 52,154,227, 99, 77, 33, 11,173, 50,147,166,127,182, 33,  1,233,221, 84, 48,123,
 198,121,237,105,155,213, 12,174,174,197,100,133,243,248, 22, 84, 12,174,206,164,
 198, 22,146,238, 91, 24,202,171,181,189,162,121, 57},
{ 85,  1, 29,145, 14,230,225, 85,194,242,140, 77,215,250,214, 40,200,226,106,  5,
 171,215,135,151, 77,226,225,111,142,246,176,156,  0,215, 18,228, 41,  8, 34,151,
  24,174,236,105, 28,  5, 39,243,194, 63,128,181, 19},
{ 44,163, 27, 35, 21, 83,238,106,156,122, 59,255,250,132, 43, 45, 12,243,  8,  9,
  16,185,194,  2,126,136,115,220,237, 47,141,167,212, 35,164, 47,217,206, 88,195,
 238, 68,125, 44,175, 49,177,138,  4,213,165,186,120},
{ 55,190, 96,216, 35,121,141,182, 26, 28,152, 34,238,248, 75,122,213,237, 99,213,
  34, 61,152,173,145,204,133,143, 64,117,119, 92,224, 76,187, 36,160,208,177, 95,
 127,213, 58,214,134, 44,121,248, 82, 63,169,191, 75},
{187, 28, 69, 29, 89,  4,160,228, 22,185, 43, 88,154, 12, 86,206, 43,199,115,152,
  40,239, 11,192, 73,228,145, 24,154, 41, 63, 49, 40, 36,224,176,100, 94, 31,100,
 152,109,111,135,185,118,207, 58, 18,247, 59,144, 33},
{117, 25, 72,154,251,194,111, 69,202,191,253,159,120,178,246, 68,171, 41,251,163,
 124,202,254,239,152, 25,  2,  5,204,223,192,231,250,120,193,179,234, 80,108,166,
 166,167,210,195, 99,135,159,118,132,143,164,128, 36},
{143, 12,156, 52,139,203,193, 61, 89,  3, 53, 84, 14,168,101,194,207, 61,113, 59,
 188, 39,200, 99, 26, 41, 88,222,211,134,178,117, 71, 15,136,150,150, 65, 88,124,
 204,128, 23, 28, 51,166,204,221,251, 63, 53, 44,190},
{203,226, 36, 10,145, 27, 54,129,243,142, 43, 63,242, 57,243, 98,229, 59, 74,201,
  41, 44, 96,199,124, 97,197, 70,118, 78,134, 66,106,138, 68,197, 64,140,187, 91,
 201, 10,138,135, 16,254,109,113,144,220,128,204, 93},
{ 29, 55,158,167,195,223,144,158,158,116, 87,219,101, 36, 71, 28,189, 52,215, 17,
 199, 92,176,139, 74,132,108,  3, 25,126, 46,191,226,239, 14,161, 44, 70,247,253,
 202,246, 58, 36, 35, 29, 77,144, 52, 14,217,139,221},
{122, 57, 40, 21, 48, 65, 99, 21, 77, 50,204, 30,233,166,117,  3, 48,  3,115,250,
 224, 78,143,108,245,144,255,199,147,114,161, 38,145, 41,107,172,132, 82, 95,202,
 166,152, 75, 83, 88,143, 25, 25,186,202,151,159,222},
{125, 19, 56,207,112, 92,184,147,239,181,113,209, 24,245,173, 57,173, 51,  3,160,
 148,255,182, 92,140,168,146,194,234, 61, 53,190,137, 15, 91,228,231,  9,111,222,
  52, 62,205,189, 90,185,129,222, 74, 19,154, 94, 29},
{161,204,117,222,253, 61,201, 66,207,106, 21,166,117,149,224,164,249, 50, 45,172,
  71,205, 29, 87,112, 81,177, 95,215,130,214,162, 83, 43,182,  9,188,112,183,111,
   5,174,231,176,103,151,117,  7,232,167, 19, 33,234},
{207,147,205, 21,140,244, 31,178,149,173,157, 33,161, 85,130,130,237,116,136, 51,
  54,137,106,123,126,234,208, 57,145, 34,116,229,209,226, 26, 86, 63,239,245,210,
  21,211, 61,189, 43, 85,215,103,160,170,234,163, 56},
{215,200,167, 19,210,166, 18, 96,224, 77,  5,145,106,148,222,103,157,196,233,132,
 109, 61,229,187,163,152, 17, 62, 27,210, 42, 67,181,  2, 23,108, 68,206,189, 76,
  58, 39,164, 43,254,  9, 87, 41, 18,228,135,212,165}
};
/* Expected decoded message (15 x 53 octets), ICD Annex C --------------------*/
static const uint8_t dec_expected[15 * 53] = {
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
/* Shift the 56-octet HAS Page out of a 64-octet C/NAV buffer at bit offset 14*/
/* (used to feed the ICD Annex C reference page through the new API).         */
static void cnav_to_has_page(const uint8_t *cnav64, uint8_t *page56)
{
    int i;
    for (i = 0; i < HAS_PAGE_BYTES; i++) {
        page56[i] = (uint8_t)getbitu(cnav64, 14 + i * 8, 8);
    }
}
/* Build a synthetic 56-octet HAS Page from a header and a 53-octet payload.  */
/* Bytes 0..2 hold the 24-bit HAS Page Header; bytes 3..55 hold the encoded   */
/* payload. The transport-level Reserved field and CRC do not exist at the    */
/* HAS layer and are not produced by this helper.                             */
static void build_page(int hass, int mt, int mid, int ms, int pid,
                       const uint8_t *enc53, uint8_t *page56)
{
    memset(page56, 0, HAS_PAGE_BYTES);
    setbitu(page56,  0, 2, hass);
    /* bits 2..3 reserved */
    setbitu(page56,  4, 2, mt);
    setbitu(page56,  6, 5, mid);
    setbitu(page56, 11, 5, ms - 1);     /* "0"=1 .. "31"=32           */
    setbitu(page56, 16, 8, pid);
    memcpy(page56 + HAS_HDR_BYTES, enc53, HAS_ENC_BYTES);
}
/* Test 1: parse the ICD Annex C page (single page, after C/NAV-to-HAS shift)*/
static int test_annexc_page(void)
{
    gal_has_t has;
    gtime_t t = {0};
    uint8_t has_page[HAS_PAGE_BYTES];
    int ret;

    has_init(&has);
    cnav_to_has_page(annexc_cnav, has_page);
    ret = has_input_page(&has, t, has_page);
    if (ret != 0) {
        fprintf(stderr, "test1: expected 0 (single page, ms=15), got %d\n",
                ret);
        return 1;
    }
    if (has.chan[0].npage != 1 || has.chan[0].mt != 1 ||
        has.chan[0].mid != 15  || has.chan[0].ms != 15 ||
        has.chan[0].pids[0] != 55) {
        fprintf(stderr, "test1: wrong header parse: npage=%d mt=%d mid=%d "
                "ms=%d pid=%d\n",
                has.chan[0].npage, has.chan[0].mt, has.chan[0].mid,
                has.chan[0].ms, has.chan[0].pids[0]);
        return 1;
    }
    /* The first byte of the encoded payload must equal the first byte of   */
    /* PID-55 in the ICD: 132. Spot-check.                                  */
    if (has.chan[0].pages[0][0] != 132) {
        fprintf(stderr, "test1: wrong encoded payload byte 0: got %u\n",
                has.chan[0].pages[0][0]);
        return 1;
    }
    printf("test1: parse Annex C HAS Page OK (header validated)\n");
    return 0;
}
/* Test 2: dummy HAS Page must be silently ignored ---------------------------*/
static int test_dummy_page(void)
{
    gal_has_t has;
    gtime_t t = {0};
    uint8_t dummy[HAS_PAGE_BYTES];
    int ret;

    has_init(&has);
    memset(dummy, 0, sizeof(dummy));
    dummy[0] = 0xAF;
    dummy[1] = 0x3B;
    dummy[2] = 0xC3;
    ret = has_input_page(&has, t, dummy);
    if (ret != 0 || has.chan[0].npage != 0) {
        fprintf(stderr, "test2: dummy page not discarded (ret=%d npage=%d)\n",
                ret, has.chan[0].npage);
        return 1;
    }
    printf("test2: dummy page discarded OK\n");
    return 0;
}
/* Test 3: end-to-end reassembly + RS decode of the 15-page ICD example ------*/
static int test_full_pipeline(void)
{
    gal_has_t has;
    gtime_t t = {0};
    uint8_t page[HAS_PAGE_BYTES];
    int order[15] = {7,3,11,0,5,12,1,9,14,4,8,2,13,6,10}; /* shuffled order */
    int ret = 0, completed = 0, i, j, mismatches = 0;

    has_init(&has);

    for (i = 0; i < 15; i++) {
        j = order[i];
        build_page(0, 1, 7, 15, pids_anc[j], enc_anc[j], page);
        ret = has_input_page(&has, t, page);
        if (ret == 1) {
            completed++;
            if (i != 14) {
                fprintf(stderr,
                    "test3: completion fired too early after %d pages\n",
                    i + 1);
                return 1;
            }
        }
    }
    if (completed != 1) {
        fprintf(stderr, "test3: completion never fired (got %d times)\n",
                completed);
        return 1;
    }
    if (has.mt != 1 || has.mid != 7 || has.ms != 15 ||
        has.n != 15 * 53) {
        fprintf(stderr, "test3: wrong assembled metadata: mt=%d mid=%d "
                "ms=%d n=%d\n", has.mt, has.mid, has.ms, has.n);
        return 1;
    }
    for (i = 0; i < 15 * 53; i++) {
        if (has.msg[i] != dec_expected[i]) {
            if (mismatches < 5) {
                fprintf(stderr, "test3: byte %d mismatch: got %u expect %u\n",
                        i, has.msg[i], dec_expected[i]);
            }
            mismatches++;
        }
    }
    if (mismatches) {
        fprintf(stderr, "test3: FAIL (%d/%d byte mismatches)\n",
                mismatches, 15 * 53);
        return 1;
    }
    printf("test3: full pipeline OK (15 shuffled HAS Pages -> %d-byte "
           "message matches ICD Annex C)\n", has.n);
    return 0;
}
int main(void)
{
    int rc = 0;
    rc |= test_annexc_page();
    rc |= test_dummy_page();
    rc |= test_full_pipeline();
    if (rc) {
        fprintf(stderr, "test_has_decode: FAIL\n");
        return 1;
    }
    printf("test_has_decode: ALL TESTS PASS\n");
    return 0;
}
