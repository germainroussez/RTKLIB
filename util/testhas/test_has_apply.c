/*------------------------------------------------------------------------------
* test_has_apply.c : standalone validation of the HAS SSR adapter
*
*          Copyright (C) 2026 by G.ROUSSEZ, All rights reserved.
*
* description :
*     Three checks are run against has_apply_corrections() in src/has.c.
*
*       (1) A complete MT1 carrying mask + orbit + clock full-set + code
*           biases + phase biases is hand-packed for two Galileo satellites
*           and two signals, written into nav.has.msg, and applied. Every
*           ssr_t field (deph, dclk, cbias, pbias, iode, t0, udi, update)
*           is checked for the expected scale-applied, sign-flipped value.
*
*       (2) MT1-slow (mask + orbit) is applied first, then MT1-fast (clock
*           only) is applied with the same Mask ID and IOD Set ID. The
*           clock corrections of the second message must reach the satel-
*           lites identified by the first one through the (Mask ID, IOD
*           Set ID) linkage cached in nav.has.
*
*       (3) Sentinel values "data not available" (HAS_NA_DR for orbit,
*           HAS_NA_DCC for clock, HAS_NA_BIAS for code/phase biases) on a
*           subset of satellites are honoured: those entries leave the
*           corresponding ssr_t fields untouched while their neighbours
*           are still applied normally.
*
* references :
*     [1] European GNSS Service Centre, Galileo High Accuracy Service
*         Signal-in-Space Interface Control Document, Issue 1.0, May 2022.
*
* usage :
*     cc -O2 -I../../src ../../src/has.c ../../src/rtkcmn.c \
*        test_has_apply.c -o test_has_apply
*     ./test_has_apply
*
* version : $Revision:$ $Date:$
* history : 2026/05/02 1.0 new
*-----------------------------------------------------------------------------*/
#include "rtklib.h"
#include <stdio.h>
#include <string.h>
#include <math.h>

/* small helper to pack a bit-stream incrementally ---------------------------*/
typedef struct {
    uint8_t buf[4096];
    int off;
} packer_t;

static void pack_u(packer_t *p, int nbits, uint32_t v)
{
    setbitu(p->buf, p->off, nbits, v);
    p->off += nbits;
}
static void pack_s(packer_t *p, int nbits, int32_t v)
{
    uint32_t u = (uint32_t)v & ((nbits == 32) ? 0xFFFFFFFFu :
                                ((1u << nbits) - 1));
    setbitu(p->buf, p->off, nbits, u);
    p->off += nbits;
}
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
/* close to almost-equal float test (in metres) ------------------------------*/
static int near_eq(double a, double b, double tol)
{
    return fabs(a - b) <= tol;
}
/* feed nav with the bit-packed MT1 currently in p, set reception time -----*/
static void install_mt1(nav_t *nav, const packer_t *p, gtime_t t)
{
    nav->has.t = t;
    nav->has.n = (p->off + 7) / 8;
    memcpy(nav->has.msg, p->buf, nav->has.n);
    /* mt/mid/ms metadata are not required by has_apply_corrections; they   */
    /* are reset here for clarity since each test fakes a fresh delivery.   */
    nav->has.mt  = 1;
    nav->has.mid = 0;
    nav->has.ms  = 0;
}
/* clear all ssr.update flags so we can detect what gets written ------------*/
static void clear_ssr(nav_t *nav)
{
    int i;
    for (i = 0; i < MAXSAT; i++) {
        memset(&nav->ssr[i], 0, sizeof(ssr_t));
    }
}
/* allocate a heap-resident nav_t (it is several MB on stack would overflow)*/
static nav_t *make_nav(void)
{
    nav_t *nav = (nav_t *)calloc(1, sizeof(nav_t));
    return nav;
}
/* pack a complete MT1 with mask + orbit + clock_full + cb + pb -------------*/
/* into p, using the test-fixed parameters below. Parameters that the       */
/* caller wants to vary (orbit values, sentinels, etc.) are passed in.      */
typedef struct {
    int toh, mask_id, iod_set_id;
    int gnss_id;                /* HAS_GNSS_GAL */
    int sat_svids[2];           /* SVIDs of the two satellites */
    int sigs[2];                /* signal indexes */
    int orbit_iod[2];
    int orbit_dr[2], orbit_dit[2], orbit_dct[2];
    int dcm;                    /* multiplier raw 0..3 -> 1..4 */
    int dcc[2];
    int cb[2][2];               /* per (sat, sig) */
    int pb[2][2];
    int pdi[2][2];
    int orbit_flag, clk_full_flag, cb_flag, pb_flag;
} mt1_fixture_t;

static void pack_mt1(packer_t *p, const mt1_fixture_t *f)
{
    int i, j;
    /* MT1 header */
    pack_u(p, 12, f->toh);
    pack_u(p, 1, 1);                            /* mask flag */
    pack_u(p, 1, f->orbit_flag);
    pack_u(p, 1, f->clk_full_flag);
    pack_u(p, 1, 0);                            /* clk sub */
    pack_u(p, 1, f->cb_flag);
    pack_u(p, 1, f->pb_flag);
    pack_u(p, 4, 0);                            /* reserved */
    pack_u(p, 5, f->mask_id);
    pack_u(p, 5, f->iod_set_id);

    /* Mask block: 1 GNSS, 2 sats, 2 sigs, no cell mask */
    pack_u(p, 4, 1);                            /* nsys */
    pack_u(p, 4, f->gnss_id);
    pack_sat_mask(p, f->sat_svids, 2);
    pack_sig_mask(p, f->sigs, 2);
    pack_u(p, 1, 0);                            /* cmaf */
    pack_u(p, 3, 0);                            /* nm */
    pack_u(p, 6, 0);                            /* mask reserved trailer */

    /* Orbit block */
    if (f->orbit_flag) {
        pack_u(p, 4, 5);                        /* VI = 5 (60 s) */
        for (i = 0; i < 2; i++) {
            pack_u(p, 10, f->orbit_iod[i]);
            pack_s(p, 13, f->orbit_dr[i]);
            pack_s(p, 12, f->orbit_dit[i]);
            pack_s(p, 12, f->orbit_dct[i]);
        }
    }
    /* Clock Full-Set */
    if (f->clk_full_flag) {
        pack_u(p, 4, 5);                        /* VI = 5 (60 s) */
        pack_u(p, 2, f->dcm);
        for (i = 0; i < 2; i++) pack_s(p, 13, f->dcc[i]);
    }
    /* Code biases (no cell mask, both cells provided) */
    if (f->cb_flag) {
        pack_u(p, 4, 4);                        /* VI = 4 (30 s) */
        for (i = 0; i < 2; i++)
            for (j = 0; j < 2; j++)
                pack_s(p, 11, f->cb[i][j]);
    }
    /* Phase biases */
    if (f->pb_flag) {
        pack_u(p, 4, 4);                        /* VI = 4 (30 s) */
        for (i = 0; i < 2; i++)
            for (j = 0; j < 2; j++) {
                pack_s(p, 11, f->pb[i][j]);
                pack_u(p,  2, f->pdi[i][j]);
            }
    }
}
/* Test 1: full MT1 round-trip ------------------------------------------------*/
static int test_full_mt1(nav_t *nav)
{
    packer_t p = {{0}, 0};
    mt1_fixture_t f;
    int sat_e1, sat_e2;
    gtime_t t0 = epoch2time((double[]){2026, 3, 6, 16, 20,  0});
    int rc, code1, code2;

    memset(&f, 0, sizeof(f));
    f.toh = 1234; f.mask_id = 7; f.iod_set_id = 11;
    f.gnss_id = HAS_GNSS_GAL;
    f.sat_svids[0] = 5;  f.sat_svids[1] = 27;
    f.sigs[0]      = 0;  f.sigs[1]      = 5;       /* E1-B and E5a-X */
    f.orbit_iod[0] = 123; f.orbit_iod[1] = 456;
    f.orbit_dr[0]  = 100; f.orbit_dr[1]  = -200;   /* +0.25 / -0.50 m */
    f.orbit_dit[0] = -50; f.orbit_dit[1] = 500;    /* -0.40 / +4.00 m */
    f.orbit_dct[0] = 10;  f.orbit_dct[1] = -10;
    f.dcm          = 1;                            /* multiplier 2     */
    f.dcc[0]       = 300; f.dcc[1] = -75;
    f.cb[0][0]     = 50;  f.cb[0][1] = -100;       /* +1.0 / -2.0 m   */
    f.cb[1][0]     = 25;  f.cb[1][1] = 0;
    f.pb[0][0]     = -30; f.pb[0][1] = 120;        /* -0.3 / +1.2 cyc */
    f.pb[1][0]     = 200; f.pb[1][1] = -200;
    f.pdi[0][0]    = 1;   f.pdi[0][1] = 2;
    f.pdi[1][0]    = 3;   f.pdi[1][1] = 0;
    f.orbit_flag = f.clk_full_flag = f.cb_flag = f.pb_flag = 1;

    pack_mt1(&p, &f);
    install_mt1(nav, &p, t0);
    has_init(&nav->has);
    nav->has.t = t0;
    nav->has.n = (p.off + 7) / 8;
    memcpy(nav->has.msg, p.buf, nav->has.n);
    clear_ssr(nav);

    rc = has_apply_corrections(nav);
    if (rc != 1) {
        fprintf(stderr, "test1: has_apply_corrections returned %d\n", rc);
        return 1;
    }
    sat_e1 = satno(SYS_GAL, 5);
    sat_e2 = satno(SYS_GAL, 27);
    if (!sat_e1 || !sat_e2) {
        fprintf(stderr, "test1: satno() failed\n"); return 1;
    }
    /* Orbit: HAS_DR=100, scale 0.0025 -> 0.25 m, sign flipped -> -0.25 */
    if (!near_eq(nav->ssr[sat_e1-1].deph[0], -0.25,  1e-9) ||
        !near_eq(nav->ssr[sat_e1-1].deph[1],  0.40,  1e-9) ||
        !near_eq(nav->ssr[sat_e1-1].deph[2], -0.08,  1e-9) ||
        !near_eq(nav->ssr[sat_e2-1].deph[0],  0.50,  1e-9) ||
        !near_eq(nav->ssr[sat_e2-1].deph[1], -4.00,  1e-9) ||
        !near_eq(nav->ssr[sat_e2-1].deph[2],  0.08,  1e-9)) {
        fprintf(stderr,
            "test1: orbit values wrong: "
            "e1=[%g %g %g] e2=[%g %g %g]\n",
            nav->ssr[sat_e1-1].deph[0], nav->ssr[sat_e1-1].deph[1],
            nav->ssr[sat_e1-1].deph[2],
            nav->ssr[sat_e2-1].deph[0], nav->ssr[sat_e2-1].deph[1],
            nav->ssr[sat_e2-1].deph[2]);
        return 1;
    }
    if (nav->ssr[sat_e1-1].iode != 123 || nav->ssr[sat_e2-1].iode != 456) {
        fprintf(stderr, "test1: iode wrong\n"); return 1;
    }
    /* Clock: dcm=2, dcc[0]=300 -> 2*300*0.0025 = 1.5 m (no sign flip) */
    if (!near_eq(nav->ssr[sat_e1-1].dclk[0],  1.50, 1e-9) ||
        !near_eq(nav->ssr[sat_e2-1].dclk[0], -0.375, 1e-9)) {
        fprintf(stderr, "test1: dclk wrong: e1=%g e2=%g\n",
                nav->ssr[sat_e1-1].dclk[0], nav->ssr[sat_e2-1].dclk[0]);
        return 1;
    }
    /* Code biases: HAS CB=50 -> 50*0.02 = 1.0 m, sign flipped -> -1.0 */
    code1 = CODE_L1B;
    code2 = CODE_L5X;
    if (!near_eq(nav->ssr[sat_e1-1].cbias[code1-1], -1.00, 1e-6) ||
        !near_eq(nav->ssr[sat_e1-1].cbias[code2-1],  2.00, 1e-6) ||
        !near_eq(nav->ssr[sat_e2-1].cbias[code1-1], -0.50, 1e-6) ||
        !near_eq(nav->ssr[sat_e2-1].cbias[code2-1],  0.00, 1e-6)) {
        fprintf(stderr,
            "test1: cbias wrong: e1[L1B]=%g e1[L5X]=%g e2[L1B]=%g e2[L5X]=%g\n",
            nav->ssr[sat_e1-1].cbias[code1-1], nav->ssr[sat_e1-1].cbias[code2-1],
            nav->ssr[sat_e2-1].cbias[code1-1], nav->ssr[sat_e2-1].cbias[code2-1]);
        return 1;
    }
    /* Phase biases: HAS PB=-30 -> -30*0.01 = -0.30 cycles, sign-flipped to */
    /* +0.30 cycles, then converted to metres via lambda = c/freq.          */
    {
        double freq_e1b = code2freq(SYS_GAL, CODE_L1B, 0);
        double freq_e5x = code2freq(SYS_GAL, CODE_L5X, 0);
        double exp_pb_e1_b = -((double)f.pb[0][0] * 0.01) * (CLIGHT / freq_e1b);
        if (!near_eq(nav->ssr[sat_e1-1].pbias[code1-1], exp_pb_e1_b, 1e-3)) {
            fprintf(stderr, "test1: pbias E1 sat e1 wrong: got=%g want=%g\n",
                    nav->ssr[sat_e1-1].pbias[code1-1], exp_pb_e1_b);
            return 1;
        }
        if (freq_e5x <= 0) {
            fprintf(stderr, "test1: code2freq(L5X) returned 0\n"); return 1;
        }
    }
    /* Update flag must be set */
    if (!nav->ssr[sat_e1-1].update || !nav->ssr[sat_e2-1].update) {
        fprintf(stderr, "test1: update flag not set\n"); return 1;
    }
    printf("test1: full MT1 round-trip OK\n");
    return 0;
}
/* Test 2: MT1-slow then MT1-fast linkage ------------------------------------*/
static int test_slow_fast_linkage(nav_t *nav)
{
    packer_t p_slow = {{0}, 0}, p_fast = {{0}, 0};
    mt1_fixture_t f;
    int sat_e1, sat_e2;
    gtime_t t0 = epoch2time((double[]){2026, 3, 6, 16, 20,  0});
    gtime_t t1 = epoch2time((double[]){2026, 3, 6, 16, 20,  5});
    int rc;

    memset(&f, 0, sizeof(f));
    f.toh = 100; f.mask_id = 12; f.iod_set_id = 19;
    f.gnss_id = HAS_GNSS_GAL;
    f.sat_svids[0] = 9;  f.sat_svids[1] = 18;
    f.sigs[0] = 0;       f.sigs[1] = 5;
    f.orbit_iod[0] = 50; f.orbit_iod[1] = 80;
    f.orbit_dr[0]  = 40; f.orbit_dr[1]  = -40;
    f.orbit_dit[0] =  0; f.orbit_dit[1] = 0;
    f.orbit_dct[0] =  0; f.orbit_dct[1] = 0;
    f.orbit_flag   = 1;

    /* slow MT1: mask + orbit only */
    pack_mt1(&p_slow, &f);

    has_init(&nav->has);
    install_mt1(nav, &p_slow, t0);
    clear_ssr(nav);
    rc = has_apply_corrections(nav);
    if (rc != 1) {
        fprintf(stderr, "test2: slow MT1 apply rc=%d\n", rc); return 1;
    }
    /* fast MT1: clock-only with same mask_id and iod_set_id */
    f.toh = 105;
    f.orbit_flag   = 0;
    f.clk_full_flag = 1;
    f.dcm  = 0;     /* multiplier 1 */
    f.dcc[0] = 200; /* 1 * 200 * 0.0025 = 0.5 m */
    f.dcc[1] = -50; /* 1 * -50 * 0.0025 = -0.125 m */

    pack_mt1(&p_fast, &f);
    install_mt1(nav, &p_fast, t1);
    rc = has_apply_corrections(nav);
    if (rc != 1) {
        fprintf(stderr, "test2: fast MT1 apply rc=%d\n", rc); return 1;
    }
    sat_e1 = satno(SYS_GAL, 9);
    sat_e2 = satno(SYS_GAL, 18);
    if (!near_eq(nav->ssr[sat_e1-1].dclk[0],  0.500, 1e-9) ||
        !near_eq(nav->ssr[sat_e2-1].dclk[0], -0.125, 1e-9)) {
        fprintf(stderr,
            "test2: linked clock wrong: sat9=%g sat18=%g\n",
            nav->ssr[sat_e1-1].dclk[0], nav->ssr[sat_e2-1].dclk[0]);
        return 1;
    }
    /* Orbit must still be present from the earlier slow MT1 */
    if (!near_eq(nav->ssr[sat_e1-1].deph[0], -0.10, 1e-9)) {
        fprintf(stderr, "test2: slow orbit lost: deph[0]=%g\n",
                nav->ssr[sat_e1-1].deph[0]);
        return 1;
    }
    printf("test2: MT1-slow + MT1-fast linkage OK\n");
    return 0;
}
/* Test 3: data-not-available sentinels are honoured ------------------------*/
static int test_na_sentinels(nav_t *nav)
{
    packer_t p = {{0}, 0};
    mt1_fixture_t f;
    int sat_a, sat_b;
    gtime_t t0 = epoch2time((double[]){2026, 3, 6, 16, 30,  0});
    int rc;

    memset(&f, 0, sizeof(f));
    f.toh = 600; f.mask_id = 4; f.iod_set_id = 5;
    f.gnss_id = HAS_GNSS_GAL;
    f.sat_svids[0] = 11; f.sat_svids[1] = 22;
    f.sigs[0] = 0;       f.sigs[1] = 5;
    f.orbit_iod[0] = 10; f.orbit_iod[1] = 20;
    /* Sat A: orbit normal, but cb / pb flagged N/A on signal 0           */
    f.orbit_dr[0]  = 80; f.orbit_dit[0] = 0; f.orbit_dct[0] = 0;
    /* Sat B: orbit N/A entirely                                         */
    f.orbit_dr[1]  = HAS_NA_DR;
    f.orbit_dit[1] = HAS_NA_DIT;
    f.orbit_dct[1] = HAS_NA_DIT;
    f.dcm = 0;
    f.dcc[0] = 100;
    f.dcc[1] = HAS_NA_DCC;
    f.cb[0][0] = HAS_NA_BIAS; f.cb[0][1] = 50;
    f.cb[1][0] = 30;          f.cb[1][1] = 30;
    f.pb[0][0] = 10;          f.pb[0][1] = HAS_NA_BIAS;
    f.pb[1][0] = 20;          f.pb[1][1] = 20;
    f.pdi[0][0] = 0; f.pdi[0][1] = 0; f.pdi[1][0] = 0; f.pdi[1][1] = 0;
    f.orbit_flag = f.clk_full_flag = f.cb_flag = f.pb_flag = 1;

    pack_mt1(&p, &f);
    has_init(&nav->has);
    install_mt1(nav, &p, t0);
    clear_ssr(nav);
    rc = has_apply_corrections(nav);
    if (rc != 1) {
        fprintf(stderr, "test3: apply rc=%d\n", rc); return 1;
    }
    sat_a = satno(SYS_GAL, 11);
    sat_b = satno(SYS_GAL, 22);

    /* Sat A orbit applied; Sat B orbit untouched (deph stays 0) */
    if (!near_eq(nav->ssr[sat_a-1].deph[0], -0.20, 1e-9)) {
        fprintf(stderr, "test3: sat A orbit wrong\n"); return 1;
    }
    if (nav->ssr[sat_b-1].deph[0] != 0.0 ||
        nav->ssr[sat_b-1].deph[1] != 0.0 ||
        nav->ssr[sat_b-1].deph[2] != 0.0) {
        fprintf(stderr,
            "test3: sat B orbit should stay zero (got %g %g %g)\n",
            nav->ssr[sat_b-1].deph[0], nav->ssr[sat_b-1].deph[1],
            nav->ssr[sat_b-1].deph[2]);
        return 1;
    }
    /* Sat A clock applied; Sat B clock untouched */
    if (!near_eq(nav->ssr[sat_a-1].dclk[0], 0.25, 1e-9)) {
        fprintf(stderr, "test3: sat A clock wrong (%g)\n",
                nav->ssr[sat_a-1].dclk[0]);
        return 1;
    }
    if (nav->ssr[sat_b-1].dclk[0] != 0.0) {
        fprintf(stderr,
            "test3: sat B clock should stay zero (%g)\n",
            nav->ssr[sat_b-1].dclk[0]);
        return 1;
    }
    /* CB on (sat_a, sig 0) is NA -> stays zero; (sat_a, sig 1) is set    */
    if (nav->ssr[sat_a-1].cbias[CODE_L1B - 1] != 0.0f) {
        fprintf(stderr, "test3: sat A L1B cbias should stay zero\n");
        return 1;
    }
    if (!near_eq(nav->ssr[sat_a-1].cbias[CODE_L5X - 1], -1.0, 1e-6)) {
        fprintf(stderr, "test3: sat A L5X cbias wrong\n"); return 1;
    }
    /* PB on (sat_a, sig 1) is NA -> stays zero */
    if (nav->ssr[sat_a-1].pbias[CODE_L5X - 1] != 0.0) {
        fprintf(stderr, "test3: sat A L5X pbias should stay zero\n");
        return 1;
    }
    printf("test3: N/A sentinels honoured OK\n");
    return 0;
}
int main(void)
{
    nav_t *nav = make_nav();
    int rc = 0;

    if (!nav) {
        fprintf(stderr, "test_has_apply: out of memory\n"); return 1;
    }
    rc |= test_full_mt1(nav);
    rc |= test_slow_fast_linkage(nav);
    rc |= test_na_sentinels(nav);

    free(nav);
    if (rc) {
        fprintf(stderr, "test_has_apply: FAIL\n");
        return 1;
    }
    printf("test_has_apply: ALL TESTS PASS\n");
    return 0;
}
