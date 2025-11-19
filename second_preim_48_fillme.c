#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>

#include "xoshiro.h"
#include <omp.h>

#define EM_TABLE_LOG2 22u
#define EM_TABLE_SIZE (1u << EM_TABLE_LOG2)
#define HT_LOG2 (EM_TABLE_LOG2 + 1u)
#define HT_SIZE (1u << HT_LOG2)

#define ROTL24_16(x) ((((x) << 16) ^ ((x) >> 8)) & 0xFFFFFF)
#define ROTL24_3(x) ((((x) << 3) ^ ((x) >> 21)) & 0xFFFFFF)

#define ROTL24_8(x) ((((x) << 8) ^ ((x) >> 16)) & 0xFFFFFF)
#define ROTL24_21(x) ((((x) << 21) ^ ((x) >> 3)) & 0xFFFFFF)

#define CV_HT_LOG2 20u
#define CV_HT_SIZE (1u << CV_HT_LOG2)

#define IV 0x010203040506ULL

static double now_sec(void)
{
	struct timespec ts;
	clock_gettime(CLOCK_MONOTONIC, &ts);
	return ts.tv_sec + ts.tv_nsec * 1e-9;
}

/*
 * the 96-bit key is stored in four 24-bit chunks in the low bits of k[0]...k[3]
 * the 48-bit plaintext is stored in two 24-bit chunks in the low bits of p[0], p[1]
 * the 48-bit ciphertext is written similarly in c
 */
void speck48_96(const uint32_t k[4], const uint32_t p[2], uint32_t c[2])
{
	uint32_t rk[23];
	uint32_t ell[3] = {k[2], k[1], k[0]};

	rk[0] = k[3];

	c[0] = p[0];
	c[1] = p[1];

	/* full key schedule */
	for (unsigned i = 0; i < 22; i++)
	{
		uint32_t new_ell = ((ROTL24_16(ell[0]) + rk[i]) ^ i) & 0xFFFFFF; // addition (+) is done mod 2**24
		rk[i + 1] = ROTL24_3(rk[i]) ^ new_ell;
		ell[0] = ell[1];
		ell[1] = ell[2];
		ell[2] = new_ell;
	}

	for (unsigned i = 0; i < 23; i++)
	{
		c[0] = (ROTL24_16(c[0]) + c[1]) & 0xFFFFFF;

		c[0] ^= rk[i];
		c[1] = ROTL24_3(c[1]) ^ c[0];
	}

	return;
}

/* the inverse cipher */
void speck48_96_inv(const uint32_t k[4], const uint32_t c[2], uint32_t p[2])
{
	uint32_t rk[23];

	uint32_t ell[3] = {k[2], k[1], k[0]};
	rk[0] = k[3];
	// full key schedule
	for (unsigned i = 0; i < 22; i++)
	{
		uint32_t new_ell = ((ROTL24_16(ell[0]) + rk[i]) ^ i) & 0xFFFFFF;
		rk[i + 1] = ROTL24_3(rk[i]) ^ new_ell;
		ell[0] = ell[1];
		ell[1] = ell[2];
		ell[2] = new_ell;
	}

	uint32_t x = c[0], y = c[1];
	for (int i = 22; i >= 0; i--)
	{
		uint32_t y0 = ROTL24_21(y ^ x);
		uint32_t x1 = (x ^ rk[i]) & 0xFFFFFF;
		uint32_t x0 = ROTL24_8((x1 - y0) & 0xFFFFFF);
		x = x0;
		y = y0;
	}

	p[0] = x & 0xFFFFFF;
	p[1] = y & 0xFFFFFF;
}

int test_sp48_inv(void)
{
	uint32_t k[4] = {0x1a1918, 0x121110, 0x0a0908, 0x020100};
	uint32_t p[2] = {0x6d2073, 0x696874};
	uint32_t c[2], d[2];

	speck48_96(k, p, c);
	printf("enc: %06X %06X\n", c[0], c[1]);
	speck48_96_inv(k, c, d);
	printf("dec: %06X %06X\n", d[0], d[1]);

	if (d[0] != p[0] || d[1] != p[1])
	{
		printf("decryption mismatch\n");
		return 0;
	}

	uint32_t ct[2] = {0x735E10, 0xB6445D};
	speck48_96_inv(k, ct, d);

	printf("dec from known ct: %06X %06X\n", d[0], d[1]);

	return 1;
}

/* Test against EP 2013/404, App. C */
bool test_vector_okay()
{
	uint32_t k[4] = {0x1a1918, 0x121110, 0x0a0908, 0x020100};
	uint32_t p[2] = {0x6d2073, 0x696874};
	uint32_t c[2];
	speck48_96(k, p, c);
	printf("%X %X\n", c[0], c[1]);

	return (c[0] == 0x735E10) && (c[1] == 0xB6445D);
}

/* The Davies-Meyer compression function based on speck48_96,
 * using an XOR feedforward
 * The input/output chaining value is given on the 48 low bits of a single 64-bit word,
 * whose 24 lower bits are set to the low half of the "plaintext"/"ciphertext" (p[0]/c[0])
 * and whose 24 higher bits are set to the high half (p[1]/c[1])
 */
uint64_t cs48_dm(const uint32_t m[4], const uint64_t h)
{
	uint32_t p[2], c[2];
	p[0] = (uint32_t)(h & 0xFFFFFF);
	p[1] = (uint32_t)((h >> 24) & 0xFFFFFF);

	speck48_96(m, p, c);

	uint32_t nh0 = (c[0] ^ p[0]) & 0xFFFFFF;

	uint32_t nh1 = (c[1] ^ p[1]) & 0xFFFFFF;
	return ((uint64_t)nh1 << 24) | (uint64_t)nh0;
}

int test_cs48_dm(void)
{
	uint32_t m[4] = {0, 1, 2, 3};
	uint64_t h = IV;
	uint64_t r = cs48_dm(m, h);
	printf("cs48_dm: %06llX\n", (unsigned long long)r);
	return r == 0x5DFD97183F91ULL;
}

/* Assumes message length is fourlen * four blocks of 24 bits, each stored as the low bits of 32-bit words
 * fourlen is stored on 48 bits (as the 48 low bits of a 64-bit word)
 * when padding is included, simply adds one block (96 bits) of padding with fourlen and zeros on higher bits
 * (The fourlen * four blocks must be considered “full”, otherwise trivial collisions are possible)
 */
uint64_t hs48(const uint32_t *m, uint64_t fourlen, int padding, int verbose)
{
	uint64_t h = IV;
	const uint32_t *mp = m;

	for (uint64_t i = 0; i < fourlen; i++)
	{
		h = cs48_dm(mp, h);
		if (verbose)
			printf("@%llu : %06X %06X %06X %06X => %06llX\n", i, mp[0], mp[1], mp[2], mp[3], h);
		mp += 4;
	}
	if (padding)
	{
		uint32_t pad[4];
		pad[0] = fourlen & 0xFFFFFF;
		pad[1] = (fourlen >> 24) & 0xFFFFFF;
		pad[2] = 0;
		pad[3] = 0;
		h = cs48_dm(pad, h);
		if (verbose)
			printf("@%llu : %06X %06X %06X %06X => %06llX\n", fourlen, pad[0], pad[1], pad[2], pad[3], h);
	}

	return h;
}

/* Computes the unique fixed-point for cs48_dm for the message m */
uint64_t get_cs48_dm_fp(uint32_t m[4])
{
	uint32_t ctz[2] = {0, 0};
	uint32_t p[2];
	speck48_96_inv(m, ctz, p);
	uint32_t h0 = p[0] & 0xFFFFFF;
	uint32_t h1 = p[1] & 0xFFFFFF;
	return ((uint64_t)h1 << 24) | (uint64_t)h0;
}

int test_cs48_dm_fp(void)
{
	uint32_t m[4] = {0, 1, 2, 3};
	uint64_t fp = get_cs48_dm_fp(m);
	uint64_t r = cs48_dm(m, fp);
	printf("fp: %06llX  F(fp): %06llX\n", (unsigned long long)fp, (unsigned long long)r);
	return r == fp;
}

static inline uint32_t rnd24(void)
{
	return (uint32_t)(xoshiro256starstar_random() & 0xFFFFFFULL);
}

static inline void random_block4(uint32_t b[4])
{
	b[0] = rnd24();
	b[1] = rnd24();
	b[2] = rnd24();
	b[3] = rnd24();
}

static inline uint64_t pack48(uint32_t lo24, uint32_t hi24)
{
	return (((uint64_t)(hi24 & 0xFFFFFF)) << 24) | (uint64_t)(lo24 & 0xFFFFFF);
}

static inline uint64_t mix64(uint64_t x)
{
	x ^= x >> 30;
	x *= 0xbf58476d1ce4e5b9ULL;
	x ^= x >> 27;
	x *= 0x94d049bb133111ebULL;
	x ^= x >> 31;
	return x;
}

struct em_entry
{
	uint64_t key_plus1;
	uint32_t m1[4];
};

static void ht_insert(struct em_entry *tab, uint64_t h, const uint32_t m1[4])
{
	uint64_t k = h + 1; // avoid 0
	uint64_t idx = mix64(h) & (HT_SIZE - 1);
	while (tab[idx].key_plus1 != 0)
	{
		idx = (idx + 1) & (HT_SIZE - 1);
	}
	tab[idx].key_plus1 = k;
	tab[idx].m1[0] = m1[0];
	tab[idx].m1[1] = m1[1];
	tab[idx].m1[2] = m1[2];
	tab[idx].m1[3] = m1[3];
}

static struct em_entry *ht_find(struct em_entry *tab, uint64_t h)
{

	uint64_t k = h + 1;
	uint64_t idx = mix64(h) & (HT_SIZE - 1);
	for (;;)
	{
		uint64_t kp1 = tab[idx].key_plus1;
		if (kp1 == 0)
			return NULL;
		if (kp1 == k)
			return &tab[idx];
		idx = (idx + 1) & (HT_SIZE - 1);
	}
}

// hellpers for attack 

struct cv_entry {
    uint64_t key_plus1;
    uint32_t idx;
};

static void cv_ht_insert(struct cv_entry *tab, uint64_t h, uint32_t idx)
{
    uint64_t k = h + 1; /* avoid 0 */
    uint64_t pos = mix64(h) & (CV_HT_SIZE - 1);

    while (tab[pos].key_plus1 != 0) {
        pos = (pos + 1) & (CV_HT_SIZE - 1);
    }
    tab[pos].key_plus1 = k;
    tab[pos].idx = idx;
}

static struct cv_entry *cv_ht_find(struct cv_entry *tab, uint64_t h)
{
    uint64_t k = h + 1;
    uint64_t pos = mix64(h) & (CV_HT_SIZE - 1);

    for (;;) {
        uint64_t kp1 = tab[pos].key_plus1;
        if (kp1 == 0)
            return NULL;
        if (kp1 == k)
            return &tab[pos];
        pos = (pos + 1) & (CV_HT_SIZE - 1);
    }
}


/* Finds a two-block expandable message for hs48, using a fixed-point
 * That is, computes m1, m2 s.t. hs48_nopad(m1||m2) = hs48_nopad(m1||m2^*),
 * where hs48_nopad is hs48 with no padding */
void find_exp_mess(uint32_t m1[4], uint32_t m2[4])
{
	static struct em_entry *T = NULL;
	if (!T)
	{
		T = (struct em_entry *)calloc(HT_SIZE, sizeof(struct em_entry));
	}

	for (uint32_t i = 0; i < EM_TABLE_SIZE; i++)
	{
		uint32_t t[4];
		random_block4(t);
		uint64_t h = cs48_dm(t, IV);
		ht_insert(T, h, t);
	}
	for (;;)
	{
		uint32_t cand[4];
		random_block4(cand);
		uint64_t fp = get_cs48_dm_fp(cand);
		struct em_entry *e = ht_find(T, fp);
		if (e)
		{
			m1[0] = e->m1[0];
			m1[1] = e->m1[1];
			m1[2] = e->m1[2];
			m1[3] = e->m1[3];
			m2[0] = cand[0];
			m2[1] = cand[1];
			m2[2] = cand[2];
			m2[3] = cand[3];
			return;
		}
	}
}

int test_em(void)
{
	double t0 = now_sec();
	uint32_t a[4], b[4];
	find_exp_mess(a, b);
	double t1 = now_sec();
	printf("find_exp_mess time: %.3f s\n", t1 - t0);

	uint64_t h1 = cs48_dm(a, IV);
	uint64_t fp = get_cs48_dm_fp(b);

	printf("EM found:\n");
	printf(" m1 = {%06X,%06X,%06X,%06X}\n", a[0], a[1], a[2], a[3]);
	printf(" m2 = {%06X,%06X,%06X,%06X}\n", b[0], b[1], b[2], b[3]);
	printf(" h  = %06llX  fp(m2) = %06llX\n",
		   (unsigned long long)h1, (unsigned long long)fp);

	if (h1 != fp)
	{
		printf("expandable message check failed: h != fp\n");
		return 0;
	}

	uint32_t buf[4 * 6]; //  m1 + up to 5 copies of m2
	for (int i = 0; i < 4; i++)
		buf[i] = a[i];

	for (int i = 0; i < 4; i++)
		buf[4 + i] = b[i]; // n=1
	uint64_t r1 = hs48(buf, 2, 0, 0);

	for (int i = 0; i < 4; i++)
		buf[8 + i] = b[i]; // n=2
	uint64_t r2 = hs48(buf, 3, 0, 0);

	for (int j = 0; j < 3; j++)
		for (int i = 0; i < 4; i++)
			buf[12 + 4 * j + i] = b[i]; // n =5
	uint64_t r5 = hs48(buf, 6, 0, 0);
	printf(" hs48(m1||m2   ) = %06llX\n", (unsigned long long)r1);
	printf(" hs48(m1||m2^2 ) = %06llX\n", (unsigned long long)r2);
	printf(" hs48(m1||m2^5 ) = %06llX\n", (unsigned long long)r5);

	return (r1 == h1) && (r2 == h1) && (r5 == h1);
}

void attack(void)
{	
	double t_attack0 = now_sec();
    const uint32_t L = 1u << 18;     //96
    const uint32_t WORDS = L * 4u;    //24
    uint32_t *mess = malloc((size_t)WORDS * sizeof(uint32_t));
    if (!mess) {
        perror("malloc mess");
        exit(1);
    }

//msg
    for (uint32_t i = 0; i < (1u << 20); i += 4) {
        mess[i + 0] = i;
        mess[i + 1] = 0;
        mess[i + 2] = 0;
        mess[i + 3] = 0;
    }

    uint64_t h_target = hs48(mess, L, 1, 0);
    printf("target hash : %06llX\n", (unsigned long long)h_target);
    if (h_target != 0x7CA651E182DBULL) {
        printf("unexpected target hash, aborting attack\n");
        free(mess);
        return;
    }


    uint64_t *chain = malloc((size_t)(L + 1) * sizeof(uint64_t));
    if (!chain) {
        perror("malloc chain");
        free(mess);
        exit(1);
    }

    chain[0] = IV;
    for (uint32_t b = 0; b < L; b++) {
        chain[b + 1] = cs48_dm(&mess[4 * b], chain[b]);
    }

    struct cv_entry *CV = calloc(CV_HT_SIZE, sizeof(struct cv_entry));
    if (!CV) {
        perror("calloc CV");
        free(chain);
        free(mess);
        exit(1);
    }

    for (uint32_t j = 0; j <= L; j++) {
        cv_ht_insert(CV, chain[j], j);
    }

    uint32_t m1[4], m2[4];
    double t0 = now_sec();
    find_exp_mess(m1, m2);
    double t1 = now_sec();
    uint64_t fp = cs48_dm(m1, IV);

    printf("Expandable message found in %.3f s\n", t1 - t0);
    printf("  m1 = {%06X,%06X,%06X,%06X}\n", m1[0], m1[1], m1[2], m1[3]);
    printf("  m2 = {%06X,%06X,%06X,%06X}\n", m2[0], m2[1], m2[2], m2[3]);
    printf("  fp = %06llX\n", (unsigned long long)fp);

    uint32_t cm[4];
    uint32_t j_match = 0;
    uint64_t trials = 0;
    double t2 = now_sec();

    for (;;) {
        random_block4(cm);
        uint64_t v = cs48_dm(cm, fp);
        struct cv_entry *e = cv_ht_find(CV, v);
        trials++;

        if (e && e->idx >= 3) { 
            j_match = e->idx;
            printf("Collision found after %llu trials, j = %u\n",
                   (unsigned long long)trials, j_match);
            printf("  v = %06llX (matches chain[%u])\n",
                   (unsigned long long)v, j_match);
            break;
        }
    }

    double t3 = now_sec();
    printf("Collision search time: %.3f s\n", t3 - t2);

    uint32_t n = j_match - 2;     
    uint32_t lenE = 1u + n;        
    printf("Using n = %u copies of m2 (len(E) = %u blocks)\n", n, lenE);

    uint32_t *mess2 = malloc((size_t)WORDS * sizeof(uint32_t));
    if (!mess2) {
        perror("malloc mess2");
        free(CV);
        free(chain);
        free(mess);
        exit(1);
    }

    size_t pos = 0;
    for (int i = 0; i < 4; i++)
        mess2[pos++] = m1[i];

    for (uint32_t k = 0; k < n; k++) {
        for (int i = 0; i < 4; i++)
            mess2[pos++] = m2[i];
    }


    for (int i = 0; i < 4; i++)

        mess2[pos++] = cm[i];

    

    for (uint32_t b = j_match; b < L; b++) {
        for (int i = 0; i < 4; i++)
            mess2[pos++] = mess[4 * b + i];
    }

    if (pos != (size_t)WORDS) {
        printf("length mismatch in constructed second preimage (pos=%zu, exp=%u)\n",
               pos, WORDS);
    }

    uint64_t h2 = hs48(mess2, L, 1, 0);
    printf("second preimage hash : %06llX\n", (unsigned long long)h2);

    if (h2 == h_target) {
        printf("attack success: hashes match\n");
    } else {
        printf("attack failed: hashes differ\n");
    }

	double t_attack1 = now_sec();
	printf("Total attack time: %.3f s\n", t_attack1 - t_attack0);
    free(mess2);
    free(CV);
    free(chain);
    free(mess);
}


int main()
{
	// random seed
	uint64_t seed[4] = {0, 1, 2, 3};
	xoshiro256starstar_random_set(seed);

	// attack();
	if (!test_vector_okay())
		return 1;
	if (!test_sp48_inv())
		return 1;
	if (!test_cs48_dm())
		return 1;
	if (!test_cs48_dm_fp())
		return 1;
	if (!test_em())
		return 1;
	
	attack();
	return 0;
}
