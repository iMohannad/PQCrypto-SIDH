#pragma comment(lib, "mpir.lib")
#include<stdio.h>
#include<math.h>
#include<gmp.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include <stdint.h>
#include <malloc.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#define BITS_TO_BYTES_CEIL(x) (((x) + 7) / 8)
#define MEM_NOT_FREE 0
#define MEM_FREE 1
#define SHAKE128_RATE 168
#define SHAKE256_RATE 136
#define NROUNDS 24
#define ROL(a, offset) ((a << offset) ^ (a >> (64-offset)))
#define PUBLIC_KEY_NPS 6
static const uint16_t sike_G = (uint16_t) 0;
static const uint16_t sike_H = (uint16_t) 1;
static const uint16_t pke_N  = (uint16_t) 2;
typedef mpz_t mp;


typedef struct {
	mp x0;
	mp x1;
} fp2;

typedef struct {
	fp2 x;
	fp2 y;
} mont_pt_t;



void mp_mod(const mp a, const mp mod, mp b) {
  mpz_mod(b, a, mod);
}

typedef unsigned char sike_msg;


typedef struct{
	// Starting curve is
	// E0/F_P2: y^2 = x^3 + x
	// ord(E0) = (2^eA * 3^eB)^2
	const char* name;

	// Prime p
	// p = lA^eA*lB^eB - 1
	const char* p;

	const char* eA;
	const char* eB;

	const char* lA;
	const char* lB;

	// Public generators for Alice: Q_A, P_A
	// Differential coordinate R_A = P_A - Q_A
	const char* xQA0;
	const char* xQA1;

	const char* yQA0;
	const char* yQA1;

	const char* xPA0;
	const char* xPA1;

	const char* yPA0;
	const char* yPA1;

	const char* xRA0;
	const char* xRA1;

	const char* yRA0;
	const char* yRA1;

	// Public generators for Bob: Q_B, P_B
	// Differential coordinate R_B = P_B - Q_B
	const char* xQB0;
	const char* xQB1;

	const char* yQB0;
	const char* yQB1;

	const char* xPB0;
	const char* xPB1;

	const char* yPB0;
	const char* yPB1;

	const char* xRB0;
	const char* xRB1;

	const char* yRB0;
	const char* yRB1;

	size_t crypto_bytes;
	size_t msg_bytes;
} sike_params_raw_t;

typedef struct _ff_Params ff_Params;

typedef struct _ff_Params {

	
	mp mod;

	void(*init)(const ff_Params *p, mp a);

	void(*add)(const ff_Params *p, const mp a, const mp b, mp c);

	void(*clear)(const ff_Params *p, mp a);

	void(*constant)(const ff_Params *p, unsigned long a, mp b);

	void(*copy)(const ff_Params *p, mp dst, const mp src);

	int(*isEqual)(const ff_Params *p, const mp a, const mp b);

	void(*invert)(const ff_Params *p, const mp a, mp b);

	int(*isBitSet)(const ff_Params *p, const mp a, const unsigned long index);

	int(*isConstant)(const ff_Params *p, const mp a, const size_t constant);

	void(*multiply)(const ff_Params *p, const mp a, const mp b, mp c);

	void(*negative)(const ff_Params *p, const mp a, mp b);

	void(*pow)(const ff_Params *p, const mp a, const mp b, mp c);

	//int(*rand)(const ff_Params *p, mp a);

	void(*square)(const ff_Params *d, const mp a, mp b);

	void(*subtract)(const ff_Params *p, const mp a, const mp b, mp c);

	void(*unity)(const ff_Params *p, mp b);

	void(*zero)(const ff_Params *p, mp a);

}tralala;



typedef struct {
	ff_Params* ffData;
	fp2 a;
	fp2 b;
	mont_pt_t P;
	mont_pt_t Q;
} mont_curve_int_t;




typedef struct {
	mont_curve_int_t EA;
	mont_curve_int_t EB;
	mp p;

	unsigned long eA;
	unsigned long lA;
	mp ordA;

	unsigned long eB;
	unsigned long lB;
	mp ordB;

	unsigned long msbA; // MSB of ordA
	unsigned long msbB; // MSB of ordB

	size_t crypto_bytes;
	size_t msg_bytes;
} sike_params_t;


static void gmp_add_fp(const ff_Params *p, const mp a, const mp b, mp c) {
  mpz_add(c, a, b);
  mpz_mod(c, c, p->mod);
}

static void gmp_constant_fp(const ff_Params *p, unsigned long a, mp b) {
  mpz_set_ui(b, a);
}

static void gmp_copy_fp(const ff_Params *p, mp dst, const mp src) {
  mpz_set(dst, src);
}

static int gmp_isequal_fp(const ff_Params *p, const mp a, const mp b) {
  return !mpz_cmp(a, b);
}

static void gmp_invert_fp(const ff_Params *p, const mp a, mp b) {
  mpz_invert(b, a, p->mod);
}

static int gmp_isBitSet_fp(const ff_Params *p, const mp a, const unsigned long index) {
  return mpz_tstbit(a, index);
}

static int gmp_isConstant(const ff_Params *p, const mp a, const size_t constant) {
  return !mpz_cmp_ui(a, constant);
}

static void gmp_multiply(const ff_Params *p, const mp a, const mp b, mp c) {
  mpz_mul(c, a, b);
  mpz_mod(c, c, p->mod);
}

static void gmp_negative(const ff_Params *p, const mp a, mp b) {
  mpz_neg(b, a);
  //gmp_printf("neg: %Zx\n", b);
  mpz_mod(b, b, p->mod);
}

static void gmp_pow(const ff_Params *p, const mp a, const mp b, mp c) {
  mpz_powm(c, a, b, p->mod);
}

size_t mp_sizeinbase(const mp a, int base) {
	return mpz_sizeinbase(a, base);
}
void mp_import(mp rop, size_t count, int order, size_t size, int endian, size_t nails, const void *op) {
  mpz_import (rop, count, order, size, endian, nails, op);
}

/*
static int lock = -1;


static __inline void delay(unsigned int count)
{
  while (count--) {}
}


int randombytes(unsigned char* random_array, unsigned long long xlen)
{ // Generation of "nbytes" of random values
  long long r, n = (long long) xlen, count = 0;

  if (lock == -1) {
    do {
      lock = open("/dev/urandom", O_RDONLY);
      if (lock == -1) {
        delay(0xFFFFF);
      }
    } while (lock == -1);
  }

  while (n > 0) {
    do {
      r = read(lock, random_array+count, (size_t) n);
      if (r == -1) {
        delay(0xFFFF);
      }
    } while (r == -1);
    count += r;
    n -= r;
  }

  return 0;
}



static int nist_rand(const ff_Params *p, mp a) {
  size_t bytes = mp_sizeinbase(p->mod, 2);
  unsigned char arr[bytes];
  randombytes(arr, bytes);
  mp_import(a, bytes, -1, 1, 1, 0, arr);
  mpz_mod(a, a, p->mod);
  return 0;
}*/

static void gmp_square(const ff_Params *d, const mp a, mp b) {
  mpz_mul(b, a, a);
  mpz_mod(b, b, d->mod);
}

static void gmp_subtract(const ff_Params *p, const mp a, const mp b, mp c) {
  mpz_sub(c, a, b);
  mpz_mod(c, c, p->mod);
}

static void gmp_unity(const ff_Params *p, mp b) {
  mpz_set_ui(b, 1);
}

static void gmp_zero(const ff_Params *p, mp a) {
  mpz_set_ui(a, 0);
}

static void gmp_init(const ff_Params *p, mp a) {
  mpz_init(a);
}

static void gmp_clear(const ff_Params* p, mp a) {
  mpz_clear(a);
}



void
set_gmp_fp_params(ff_Params *params) {
	params->init = gmp_init;
	params->add = gmp_add_fp;
	params->clear = gmp_clear;
	params->constant = gmp_constant_fp;
	params->copy = gmp_copy_fp;
	params->isEqual = gmp_isequal_fp;
	params->invert = gmp_invert_fp;
	params->isBitSet = gmp_isBitSet_fp;
	params->isConstant = gmp_isConstant;
	params->multiply = gmp_multiply;
	params->negative = gmp_negative;
	params->pow = gmp_pow;
	//params->rand = nist_rand;
	params->square = gmp_square;
	params->subtract = gmp_subtract;
	params->unity = gmp_unity;
	params->zero = gmp_zero;
};

void fp_Init(const ff_Params* p, mp a) {
  p->init(p, a);
}

void fp_Add(const ff_Params *p, const mp a, const mp b, mp c) {
  p->add(p, a, b, c);
}

void fp_Clear(const ff_Params *p, mp a) {
  p->clear(p, a);
}

void fp_Constant(const ff_Params *p, unsigned long a, mp b) {
  p->constant(p, a, b);
}

void fp_Copy(const ff_Params *p, mp dst, const mp src) {
  p->copy(p, dst, src);
}

int fp_IsEqual(const ff_Params *p, const mp a, const mp b) {
  return p->isEqual(p, a, b);
}

void fp_Invert(const ff_Params *p, const mp a, mp b) {
  p->invert(p, a, b);
}

int fp_IsBitSet(const ff_Params *p, const mp a, const unsigned long index) {
  return p->isBitSet(p, a, index);
}

int fp_IsConstant(const ff_Params *p, const mp a, const size_t constant) {
  return p->isConstant(p, a, constant);
}

void fp_Multiply(const ff_Params *p, const mp a, const mp b, mp c) {
  p->multiply(p, a, b, c);
}

void fp_Negative(const ff_Params *p, const mp a, mp b) {
  p->negative(p, a, b);
}

void fp_Pow(const ff_Params *p, const mp a, const mp b, mp c) {
  p->pow(p, a, b, c);
}

/*void fp_Rand(const ff_Params *p, mp a) {
  p->rand(p, a);
}*/


void fp_Square(const ff_Params *p, const mp a, mp b) {
  p->square(p, a, b);
}

void fp_Subtract(const ff_Params *p, const mp a, const mp b, mp c) {
  p->subtract(p, a, b, c);
}

void fp_Unity(const ff_Params *p, mp b) {
  p->unity(p, b);
}

void fp_Zero(const ff_Params *p, mp a) {
  p->zero(p, a);
}

void fp_ImportHex(const char *hexStr, mp a) {
  mpz_set_str(a, hexStr, 0);
}

void mp_pow(const unsigned long a, const unsigned long b, mp c) {
	mpz_ui_pow_ui(c, a, b);
}


void
fp2_Init(const ff_Params* p, fp2* fp2)
{
	fp_Init(p, fp2->x0);
	fp_Init(p, fp2->x1);
}

void
fp2_Init_set(const ff_Params* p, fp2* fp2, unsigned long x0, unsigned long x1)
{
	fp2_Init(p, fp2);
	fp_Constant(p, x0, fp2->x0);
	fp_Constant(p, x1, fp2->x1);
}

void sike_setup_params(const sike_params_raw_t *raw, sike_params_t *params) {
	// Base curve -> Coefficients are null
	mont_curve_int_t* EA = &params->EA;
	mont_curve_int_t* EB = &params->EB;

	params->eA = (unsigned long)strtol(raw->eA, NULL, 0);
	params->eB = (unsigned long)strtol(raw->eB, NULL, 0);

	params->lA = (unsigned long)strtol(raw->lA, NULL, 0);
	params->lB = (unsigned long)strtol(raw->lB, NULL, 0);

	ff_Params* ffpA = (ff_Params*)malloc(sizeof(ff_Params));
	ff_Params* ffpB = (ff_Params*)malloc(sizeof(ff_Params));


	set_gmp_fp_params(ffpA);
	set_gmp_fp_params(ffpB);

	EA->ffData = ffpA;

	fp_Init(ffpA, params->p);
	fp_ImportHex(raw->p, params->p);

	fp_Init(ffpA, params->ordA);

	mp_pow(params->lA, params->eA, params->ordA);
	params->msbA = mp_sizeinbase(params->ordA, 2);

	fp_Init(ffpA, EA->ffData->mod);

	fp_Init(ffpA, EA->P.x.x0);
	fp_Init(ffpA, EA->P.x.x1);
	fp_Init(ffpA, EA->P.y.x0);
	fp_Init(ffpA, EA->P.y.x1);
	fp_Init(ffpA, EA->Q.x.x0);
	fp_Init(ffpA, EA->Q.x.x1);
	fp_Init(ffpA, EA->Q.y.x0);
	fp_Init(ffpA, EA->Q.y.x1);

	fp_ImportHex(raw->p, EA->ffData->mod);


	fp_ImportHex(raw->xPA0, EA->P.x.x0);
	fp_ImportHex(raw->xPA1, EA->P.x.x1);
	fp_ImportHex(raw->yPA0, EA->P.y.x0);
	fp_ImportHex(raw->yPA1, EA->P.y.x1);
	fp_ImportHex(raw->xQA0, EA->Q.x.x0);
	fp_ImportHex(raw->xQA1, EA->Q.x.x1);
	fp_ImportHex(raw->yQA0, EA->Q.y.x0);
	fp_ImportHex(raw->yQA1, EA->Q.y.x1);

	fp2_Init_set(ffpA, &EA->a, 0, 0);
	fp2_Init_set(ffpB, &EA->b, 1, 0);

	EB->ffData = ffpB;

	fp_Init(ffpB, params->ordB);

	mp_pow(params->lB, params->eB, params->ordB);
	params->msbB = mp_sizeinbase(params->ordB, 2);

	fp_Init(ffpB, EB->ffData->mod);

	fp_Init(ffpB, EB->P.x.x0);
	fp_Init(ffpB, EB->P.x.x1);
	fp_Init(ffpB, EB->P.y.x0);
	fp_Init(ffpB, EB->P.y.x1);
	fp_Init(ffpB, EB->Q.x.x0);
	fp_Init(ffpB, EB->Q.x.x1);
	fp_Init(ffpB, EB->Q.y.x0);
	fp_Init(ffpB, EB->Q.y.x1);

	fp_ImportHex(raw->p, EB->ffData->mod);

	fp_ImportHex(raw->xPB0, EB->P.x.x0);
	fp_ImportHex(raw->xPB1, EB->P.x.x1);
	fp_ImportHex(raw->yPB0, EB->P.y.x0);
	fp_ImportHex(raw->yPB1, EB->P.y.x1);
	fp_ImportHex(raw->xQB0, EB->Q.x.x0);
	fp_ImportHex(raw->xQB1, EB->Q.x.x1);
	fp_ImportHex(raw->yQB0, EB->Q.y.x0);
	fp_ImportHex(raw->yQB1, EB->Q.y.x1);

	fp2_Init_set(ffpB, &EB->a, 0, 0);
	fp2_Init_set(ffpB, &EB->b, 1, 0);

	params->crypto_bytes = raw->crypto_bytes;
	params->msg_bytes = raw->msg_bytes;
}

sike_params_raw_t SIKEp751 = {
	.name = "SIKEp751",
	.p = "0x6FE5D541F71C0E12909F97BADC668562B5045CB25748084E9867D6EBE876DA959B1A13F7CC76E3EC968549F878A8EEAFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF",
	.eA = "372",
	.eB = "239",
	.lA = "2",
	.lB = "3",
	.xQA0 = "0x3E82027A38E9429C8D36FF46BCC93FA23F89F6BE06D2B1317AD90438621783FDB7A4AD3E83E86CAE096D5DB822C98E561E008FA0E3F3B9AC2F40C56D6FA4A58A20449AF1F1335661D14AB7347693632646086CE3ACD54B0346F5CCE233E9",
	.xQA1 = "0x0",
	.yQA0 = "0x3BBF8DCD4E7EB6236F5F598D56EB5E15915A755883B7C331B043DA010E6A163A7421DFA8378D1E911F50BF3F721A8ED5950D80325A8D0F147EF3BD0CFEC5236C7FAC9E69F7FD5A99EBEC3B5B8B000F8EEA737089343012E0D620BFB341D5",
	.yQA1 = "0x0",
	.xPA0 = "0x54921C31F0DC9531CB890FC5EC66DF2E7F0D55761363C6E375DA69B0682CABE5C0FFFCBE6E1AD46563F042FA06B9F207FCF3CDD2673652828FF50C3F7B755C0BE072950D16CA747C146775C0267A401FFC738B03A49E9A36B39572AFB363",
	.xPA1 = "0x28849BC0D81E01993137A5B63D6E633C4E97AB4FF118CCF63DFE623092AC86B6D4A9B751797CBA1A177500E9EB5AF7852B7DF02C334844D652EFC4729178A1DBAD8CA47BB7E757C6D43B799811A63BEBE649C18101F03AD752CDCD73BF66",
	.yPA0 = "0x196119D87272DC3AA7223476C8C3269D48CAEFAE692F68DCF2D6E1BEB5B97525D5026C157C7C740B41ADE80A8CF2E1E0B37E5F5FD4ED88235BF7404BE39189C137E21C035EF6339D7FACBA38E72D69043710E76266A5FC14EFB95E5FBC7C",
	.yPA1 = "0xD3AC09A67D59CC8D78B0FA6681AE78BDF0C8F558E3866005E4355B0B199318D9CDD67C0A7DB234F9EA1EC4C5F1E59168B7DBD14281F09E8DF904A3D574CAD526DC5A3667490ADE1A4C13B09F7B115C4E488FD4DD5F7670B5897322AD41DA",
	.xRA0 = "0x22A0B5A35A2B0C56135A7CEC5CFB97964A7C6226FE909F374362A8ECA3AB14A1B7B0C87AC875DCE5888D83B623BF0011A4AC138F62EF6B2D2D84F636548A9F920F238336E5A36E45E4055940E3C94385B8FC5374396432EEF2AE178CEFDB",
	.xRA1 = "0xF9C4AFCDA809C3358B096B250C69B20310FDF2EF631711AA4EFEC49A4E76483F320B793F2EBC63365EED14AA3F6EA33FEB56796F011BA6C6DFB4D0A00AAC4D2786646D914AD026CBB4A592EC74B5485372E51382D44528DD491B83D95471",
        .yRA0 =
"0x77B3BB69009428A327D43CA60169715F547454F88CD017B32DF58A7252C2B3C3D00D52CCD3133D54041D8BCAEA291F2057202328712CD395575CD7CCD3CE70C0A1EBF633BA946559458878F41F9FDD1727E2C31125B2FE5B713067048299",
        .yRA1 =
"0x28849BC0D81E01993137A5B63D6E633C4E97AB4FF118CCF63DFE623092AC86B6D4A9B751797CBA1A177500E9EB5AF7852B7DF02C334844D652EFC4729178A1DBAD8CA47BB7E757C6D43B799811A63BEBE649C18101F03AD752CDCD73BF66",
	.xQB0 = "0x2F1D80EF06EF960A01AB8FF409A2F8D5BCE859ED725DE145FE2D525160E0A3AD8E17B9F9238CD5E69CF26DF237429BD3778659023B9ECB610E30288A7770D3785AAAA4D646C576AECB94B919AEEDD9E1DF566C1D26D376ED2325DCC93103",
	.xQB1 = "0x77B3BB69009428A327D43CA60169715F547454F88CD213452DF58A7252C2B3C3D00D52CCD3133D54041D8BCAEA291F2057202328712CD395575CD7CCD3CE70C0A1EBF633BA946559458878F41F9FDD1727E2C31125B2FE5B713067093331",
	.yQB0 =  "0x127A46D082A1ACAF351F09AB55A15445287ED1CC55DC35892123951D4B6E302C5129C049EEB399A6EDB2EEB2F9B0A94F06CDFB3EADE76EBA0C8419745E97D12754F00E898A315B529122CFE3CA6BBC6BAF5F6BA40BB91479226A06878949",
	.yQB1 = 
"0x00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000",
	.xPB0 = "0x5FD1A3C4DD0F630974196FED3519152BC7098B9E2B121ECA46BD10A5CC9F4BCC6C689B8E4C063B3798075FCEE6EDAA9EB108B3CD00495CF04DD8CE4A08FBE685A127D40E45F4CF45098A578DEB44368699394C43BFC9BC5E00052F78E8DF",
	.xPB1 = "0x2B88A03360B3389547732C9140C05DEA6516881FE108211BE887CC43FCB80C06A1D86FF5457D3BB7DB936394EC33821AA39333A60AF84B537974CFA0BA8287D699D2BF79BA559026C64A6ED610501D2357C10B9A6C8F837424922275ACBF",
	.yPB0 = "0x53B55053E3F04FC315EFB1B7B2C4AFCB4FEF12CE744AF3B243C6E6B1417E94A78D4980DDE181896464923E01AACC3DA040A0747CA67554A352684DA207C49022D930732DF6BD0BF37E1F5C16917669A70F88059C1C739A79D7CFA0C529D9",
	.yPB1 = "0x44E44196909252ECD7B9164323815294F02AED22C4E4EB43D2CE2BC5F29EB575D45CA8B6B4C4242E369AE3A1EFC844E9D1C57B0AE3374BC2CEDAD16B0C699158332E2D9AB3F0025C0348C5F70FDC4DD7C4865E64B8B843F03D807447D5E2",
	.xRB0 = "0x77B3BB69009428A327D43CA60169715F547454F88CD017B32DF58A7252C2B3C3D00D52CCD3133D54041D8BCAEA291F2057202328712CD395575CD7CCD3CE70C0A1EBF633BA946559458878F41F9FDD1727E2C31125B2FE5B713067048297",
	.xRB1 = "0x6D91393A57DBF47FD6DCF841F17ECD719CAE1D33C6832A75B0F168855BCC38D2A4792DFF9BC86DEACA10B1AA808D539B167D73BBA32168687FA3F85AE93A1ADDE5BD1FD5B681DCC6C34454D4496976C22D80C95E42B12576FC0FB4074B9F",
	.yRB0 = "0x53B55053E3F04FC315EFB1B7B2C4AFCB4FEF12CE744AF3B243C6E6B1417E94A78D4980DDE181896464923E01AACC3DA040A0747CA67554A352684DA207C49022D930732DF6BD0BF37E1F5C16917669A70F88059C1C739A79D7CFA0C529D9",
	.yRB1 = "0x44E44196909252ECD7B9164323815294F02AED22C4E4EB43D2CE2BC5F29EB575D45CA8B6B4C4242E369AE3A1EFC844E9D1C57B0AE3374BC2CEDAD16B0C699158332E2D9AB3F0025C0348C5F70FDC4DD7C4865E64B8B843F03D807447D5EF",
	.crypto_bytes = 24,
	.msg_bytes = 32
 };

 typedef mp sike_private_key;
 typedef struct {
  ff_Params* ffData;
  fp2 xP;
  fp2 xQ;
  fp2 xR;
} sike_public_key_t;

typedef enum {ALICE, BOB} party_t;

void public_key_init(ff_Params* p, sike_public_key_t* pk) {
  pk->ffData = p;
  fp2_Init(p, &pk->xP);
  fp2_Init(p, &pk->xQ);
  fp2_Init(p, &pk->xR);
}

void ostoi(const unsigned char* to_dec, size_t to_decLen, mp dec) {
  mp_import(dec, to_decLen, -1, 1, 1, 0, to_dec);
}

void clear_free( void* ptr, size_t size, int free_mem )
{
  if (ptr) {
    memset(ptr, 0, size);
    *(volatile unsigned char*)ptr = *(volatile unsigned char*)ptr;
    if (free_mem == MEM_FREE)
      free(ptr);
  }
}

void
fp2_Copy( const ff_Params* p, const fp2* a, fp2* b )
{
  fp_Copy( p, b->x0, a->x0 );
  fp_Copy( p, b->x1, a->x1 );
}



int
fp2_IsEqual( const ff_Params* p, const fp2* a1, const fp2* a2 )
{
  return fp_IsEqual(p, a1->x0, a2->x0) && fp_IsEqual(p, a1->x1, a2->x1);
}



void
fp2_Set( const ff_Params* p, fp2* fp2, unsigned long x0, unsigned long x1 )
{
  fp_Constant(p, x0, fp2->x0);
  fp_Constant(p, x1, fp2->x1);
}


void
fp2_Invert( const ff_Params* p, const fp2* a, fp2* b )
{
  mp mul0;
  mp mul1;
  fp_Init(p, mul0);
  fp_Init(p, mul1);

  fp_Square( p, a->x0, mul0 );
  fp_Square( p, a->x1, mul1 );

  fp_Add( p, mul0, mul1, mul0 );
  fp_Invert( p, mul0, mul0 );

  fp_Negative( p, a->x1, mul1 );

  fp_Multiply( p, a->x0, mul0, b->x0 );
  fp_Multiply( p, mul1, mul0, b->x1 );

  fp_Clear(p, mul0);
  fp_Clear(p, mul1);
}



void
fp2_Multiply( const ff_Params* p,
              const fp2*  a,
              const fp2*  b,
              fp2*  c )
{
  mp mul0;
  mp mul1;
  mp adda;
  mp addb;

  fp_Init(p, mul0);
  fp_Init(p, mul1);
  fp_Init(p, adda);
  fp_Init(p, addb);

  fp_Multiply( p, a->x0, b->x0, mul0 );
  fp_Multiply( p, a->x1, b->x1, mul1 );

  fp_Add( p, a->x0, a->x1, adda );
  fp_Add( p, b->x0, b->x1, addb );

  fp_Subtract( p, mul0, mul1, c->x0 );

  fp_Add( p, mul0, mul1, mul0 );
  fp_Multiply( p, adda, addb, mul1 );

  fp_Subtract( p, mul1, mul0, c->x1 );

  fp_Clear(p, mul0);
  fp_Clear(p, mul1);
  fp_Clear(p, adda);
  fp_Clear(p, addb);
}


void
fp2_Negative( const ff_Params* p, const fp2* a, fp2* b )
{
  fp_Negative( p, a->x0, b->x0 );
  fp_Negative( p, a->x1, b->x1 );
}


void
fp2_Square( const ff_Params* p, const fp2* a, fp2* b )
{

  fp2_Multiply( p, a, a, b );
}


void
fp2_Sub( const ff_Params* p, const fp2* a, const fp2* b, fp2* c )
{
  fp_Subtract( p, a->x0, b->x0, c->x0 );
  fp_Subtract( p, a->x1, b->x1, c->x1 );
}


void
fp2_AddDeg1( const ff_Params* p, const fp2* a, const fp2* b, fp2* c)
{
  fp_Add( p, a->x0, b->x0, c->x0 );
}

void
fp2_Add( const ff_Params* p, const fp2* a, const fp2* b, fp2* c )
{
  fp_Add(p, a->x0, b->x0, c->x0);
  fp_Add(p, a->x1, b->x1, c->x1);
}
void fp2_Clear( const ff_Params* p, fp2* fp2)
{
  fp_Clear(p, fp2->x0);
  fp_Clear(p, fp2->x1);
}

/*void
fp2_Rand( const ff_Params* p, fp2* a )
{
  fp_Rand(p, a->x0);
  fp_Rand(p, a->x1);
}
*/
int
fp2_IsConst( const ff_Params* p, const fp2* a, unsigned long x0, unsigned long x1 ) {
  return p->isConstant(p, a->x0, x0) && p->isConstant(p, a->x1, x1);
}


void
fp2_Sqrt( const ff_Params* p, const fp2* a, fp2* b, int sol)
{
  mp t0, t1, t2, t3, p14, p34, inv2;
  fp_Init(p, t0);
  fp_Init(p, t1);
  fp_Init(p, t2);
  fp_Init(p, t3);
  fp_Init(p, p14);
  fp_Init(p, p34);
  fp_Init(p, inv2);

  fp_Constant(p, 2, inv2);

  // (p + 1) / 4
  mpz_add_ui(p14, p->mod, 1);
  mpz_tdiv_q_2exp(p14, p14, 2);

  // (p - 3) / 4
  mpz_sub_ui(p34, p->mod, 3);
  mpz_tdiv_q_2exp(p34, p34, 2);

  fp_Invert(p, inv2, inv2);
  fp_Square(p, a->x0, t0);
  fp_Square(p, a->x1, t1);
  fp_Add(p, t0, t1, t0);
  fp_Pow(p, t0, p14, t1);
  fp_Add(p, a->x0, t1, t0);
  fp_Multiply(p, t0, inv2, t0);
  //p->half(p, t0);
  fp_Pow(p, t0, p34, t2);

  //fp_Multiply(p, t0, t2, t1);
  fp_Pow(p, t0, p14, t1);

  fp_Multiply(p, t2, a->x1, t2);
  fp_Multiply(p, t2, inv2, t2);
  fp_Square(p, t1, t3);

  if (fp_IsEqual(p, t3, t0)) {
    if (sol == 0) {
      fp_Copy(p, b->x0, t1);
      fp_Copy(p, b->x1, t2);
    } else if (sol == 1) {
      fp_Negative(p, t1, b->x0);
      fp_Negative(p, t2, b->x1);
    }
  } else {
    if (sol == 0) {
      fp_Copy(p, b->x0, t2);
      fp_Negative(p, t1, b->x1);
    } else {
      fp_Negative(p, t2, b->x0);
      fp_Copy(p, b->x1, t1);
    }
  }

  fp_Clear(p, t0);
  fp_Clear(p, t1);
  fp_Clear(p, t2);
  fp_Clear(p, t3);
  fp_Clear(p, p14);
  fp_Clear(p, p34);
  fp_Clear(p, inv2);
}

static void mont_set_inf_affine(const mont_curve_int_t* curve, mont_pt_t *P) {

  const ff_Params *p = curve->ffData;
  fp2_Set( p, &P->x, 0, 0 );
  fp2_Set( p, &P->y, 0, 0 );
}
static int mont_is_inf_affine(const mont_curve_int_t* curve, const mont_pt_t *P) {
  return fp2_IsConst( curve->ffData, &P->y, 0, 0 );
}


void xDBL(const mont_curve_int_t *curve, const mont_pt_t *P, mont_pt_t *R) {
  const ff_Params *p = curve->ffData;

  const fp2* a = &curve->a;
  const fp2* b = &curve->b;

  // x3 = b*(3*x1^2+2*a*x1+1)^2/(2*b*y1)^2-a-x1-x1
  // y3 = (2*x1+x1+a)*(3*x1^2+2*a*x1+1)/(2*b*y1)-b*(3*x1^2+2*a*x1+1)^3/(2*b*y1)^3-y1

  if (mont_is_inf_affine(curve, P) ) {
    mont_set_inf_affine(curve, R);
    return;
  }

  fp2 t0 = { 0 }, t1 = { 0 }, t2 = { 0 };
  fp2_Init(p, &t0);
  fp2_Init(p, &t1);
  fp2_Init(p, &t2);

  fp2_Set(p, &t2, 1, 0);                        // t2 = 1

  fp2_Square(p, &P->x, &t0);                    // t0 = x1^2
  fp2_Add(p, &t0, &t0, &t1);                    // t0 = 2*x1^2
  fp2_Add(p, &t0, &t1, &t0);                    // t0 = 3*x1^2

  fp2_Multiply(p, a, &P->x, &t1);               // t1 = a*x1
  fp2_Add(p, &t1, &t1, &t1);                    // t1 = 2*a*x1


  fp2_Add(p, &t0, &t1, &t0);                    // t0 = 3*x1^2+2*a*x1
  fp2_Add(p, &t0, &t2, &t0);                    // t0 = 3*x1^2+2*a*x1+1

  fp2_Multiply(p, b, &P->y, &t1);               // t1 = b*y1
  fp2_Add(p, &t1, &t1, &t1);                    // t1 = 2*b*y1
  fp2_Invert(p, &t1, &t1);                      // t1 = 1 / (2*b*y1)

  fp2_Multiply(p, &t0, &t1, &t0);               // t0 = (3*x1^2+2*a*x1+1) / (2*b*y1)

  fp2_Square(p, &t0, &t1);                      // t1 = (3*x1^2+2*a*x1+1)^2 / (2*b*y1)^2

  fp2_Multiply(p, b, &t1, &t2);                 // t2 = b*(3*x1^2+2*a*x1+1)^2 / (2*b*y1)^2
  fp2_Sub(p, &t2, a, &t2);                      // t2 = b*(3*x1^2+2*a*x1+1)^2 / (2*b*y1)^2 - a
  fp2_Sub(p, &t2, &P->x, &t2);                  // t2 = b*(3*x1^2+2*a*x1+1)^2 / (2*b*y1)^2 - a - x1
  fp2_Sub(p, &t2, &P->x, &t2);                  // t2 = b*(3*x1^2+2*a*x1+1)^2 / (2*b*y1)^2 - a - x1 - x1

  fp2_Multiply(p, &t0, &t1, &t1);               // t1 = (3*x1^2+2*a*x1+1)^3 / (2*b*y1)^3
  fp2_Multiply(p, b, &t1, &t1);                 // t1 = b*(3*x1^2+2*a*x1+1)^3 / (2*b*y1)^3
  fp2_Add(p, &t1, &P->y, &t1);                  // t1 = b*(3*x1^2+2*a*x1+1)^3 / (2*b*y1)^3 + y1

  fp2_Add(p, &P->x, &P->x, &R->y);              // x3 = 2*x1
  fp2_Add(p, &R->y, &P->x, &R->y);              // y3 = 2*x1+x1
  fp2_Add(p, &R->y, a, &R->y);                  // y3 = 2*x1+x1+a
  fp2_Multiply(p, &R->y, &t0, &R->y);           // y3 = (2*x1+x1+a)*(3*x1^2+2*a*x1+1)/(2*b*y1)
  fp2_Sub(p, &R->y, &t1, &R->y);                // y3 = (2*x1+x1+a)*(3*x1^2+2*a*x1+1)/(2*b*y1) - (b*(3*x1^2+2*a*x1+1)^3/(2*b*y1)^3 + y1)

  fp2_Copy(p, &t2, &R->x);                      // x3 = b*(3*x1^2+2*a*x1+1)^2 / (2*b*y1)^2 - a - x1 - x1

  fp2_Clear(p, &t0);
  fp2_Clear(p, &t1);
  fp2_Clear(p, &t2);

}

void mont_pt_copy(const ff_Params* p, const mont_pt_t* src, mont_pt_t* dst) {
  if (src != dst) {
    fp2_Copy(p, &src->x, &dst->x);
    fp2_Copy(p, &src->y, &dst->y);
  }
}


void xDBLe(const mont_curve_int_t *curve, const mont_pt_t *P, int e, mont_pt_t *R) {
  mont_pt_copy(curve->ffData, P, R);
  /*printf("Values that are copied now P = \n");
 gmp_printf("%Zx \n",(*P).x.x0);
  printf("Values that are copied now R = \n");
 gmp_printf("%Zx \n",(*R).x.x0);
  printf("The properties of the curve are -- \n");
  printf("\n");
  printf("P (mont_pt_t)= \n");
  gmp_printf(" %Zx \n",(curve->P).x.x0);
  printf("\n");
  printf("Q (mont_pt_t)= \n");
  gmp_printf(" %Zx \n",(curve->Q).x.x0);
printf("\n");
  printf("a (fp2)= \n");
  gmp_printf(" %Zx \n",(curve->a).x0);
printf("\n");
  printf("b (fp2)= \n");
  gmp_printf(" %Zx \n",(curve->b).x0);
 */
  for (int j = 0; j < e; ++j)
    xDBL(curve, R, R);
}

void xADD(const mont_curve_int_t *curve, const mont_pt_t *P, const mont_pt_t *Q, mont_pt_t *R) {

  // x3 = b*(y2-y1)^2/(x2-x1)^2-a-x1-x2
  // y3 = (2*x1+x2+a)*(y2-y1)/(x2-x1)-b*(y2-y1)^3/(x2-x1)^3-y1
  // y3 = ((2*x1)+x2+a) * ((y2-y1)/(x2-x1)) - b*((y2-y1)^3/(x2-x1)^3) - y1

  const ff_Params *p = curve->ffData;

  const fp2* a = &curve->a;
  const fp2* b = &curve->b;

  fp2 t0 = { 0 }, t1 = { 0 }, t2 = { 0 };
  fp2_Init(p, &t0);
  fp2_Init(p, &t1);
  fp2_Init(p, &t2);

  fp2_Negative(p, &Q->y, &t0);

  if (mont_is_inf_affine(curve, P) ) {
    mont_pt_copy(p, Q, R);
  }
  else if (mont_is_inf_affine(curve, Q) ) {
    mont_pt_copy(p, P, R);
  }
  else if ( fp2_IsEqual( p, &P->x, &Q->x) && fp2_IsEqual( p, &P->y, &Q->y ) ) {
   
    xDBL(curve, P, R);
  }
  else if ( fp2_IsEqual( p, &P->x, &Q->x) && fp2_IsEqual(p, &P->y, &t0) ) {
    
    mont_set_inf_affine(curve, R);
  }
  else {
    

    fp2_Sub(p, &Q->y, &P->y, &t0);               // t0 = y2-y1
    fp2_Sub(p, &Q->x, &P->x, &t1);               // t1 = x2-x1
    fp2_Invert(p, &t1, &t1);                     // t1 = 1/(x2-x1)
    fp2_Multiply(p, &t0, &t1, &t0);              // t0 = (y2-y1)/(x2-x1)

    fp2_Square(p, &t0, &t1);                     // t1 = (y2-y1)^2/(x2-x1)^2

    fp2_Add(p, &P->x, &P->x, &t2);               // t2 = 2*x1
    fp2_Add(p, &t2, &Q->x, &t2);                 // t2 = 2*x1+x2
    fp2_Add(p, &t2, a, &t2);                     // t2 = 2*x1+x2+a
    fp2_Multiply(p, &t2, &t0, &t2);              // t2 = (2*x1+x2+a)*(y2-y1)/(x2-x1)

    fp2_Multiply(p, &t0, &t1, &t0);              // t0 = (y2-y1)^3/(x2-x1)^3
    fp2_Multiply(p, b, &t0, &t0);                // t0 = b*(y2-y1)^3/(x2-x1)^3
    fp2_Add(p, &t0, &P->y, &t0);                 // t0 = b*(y2-y1)^3/(x2-x1)^3+y1

    fp2_Sub(p, &t2, &t0, &t0);                   // t0 = (2*x1+x2+a)*(y2-y1)/(x2-x1)-b*(y2-y1)^3/(x2-x1)^3-y1

    fp2_Multiply(p, b, &t1, &t1);                // t1 = b*(y2-y1)^2/(x2-x1)^2
    fp2_Sub(p, &t1, a, &t1);                     // t1 = b*(y2-y1)^2/(x2-x1)^2-a
    fp2_Sub(p, &t1, &P->x, &t1);                 // t1 = b*(y2-y1)^2/(x2-x1)^2-a-x1

    fp2_Sub(p, &t1, &Q->x, &R->x);               // x3 = b*(y2-y1)^2/(x2-x1)^2-a-x1-x2

    fp2_Copy(p, &t0, &R->y);                     // y3 = (2*x1+x2+a)*(y2-y1)/(x2-x1)-(b*(y2-y1)^3/(x2-x1)^3+y1)

  }

  fp2_Clear(p, &t0);
  fp2_Clear(p, &t1);
  fp2_Clear(p, &t2);
}

void mont_pt_init(const ff_Params* p, mont_pt_t* pt) {
  fp2_Init(p, &pt->x);
  fp2_Init(p, &pt->y);
}

void mont_pt_clear(const ff_Params* p, mont_pt_t* pt) {
  fp2_Clear(p, &pt->x);
  fp2_Clear(p, &pt->y);
}

void xTPL(const mont_curve_int_t *curve, const mont_pt_t *P, mont_pt_t *R) {
  const ff_Params* p = curve->ffData;
  mont_pt_t T = { 0 };
  mont_pt_init(p, &T);

  xDBL(curve, P, &T);
  xADD(curve, P, &T, R);

  mont_pt_clear(p, &T);
}


void xTPLe(const mont_curve_int_t *curve, const mont_pt_t *P, int e, mont_pt_t *R) {
  mont_pt_copy(curve->ffData, P, R);
  for (int j = 0; j < e; ++j)
    xTPL(curve, R, R);
}

void mont_curve_copy(const ff_Params* p, const mont_curve_int_t* curve, mont_curve_int_t* curvecopy) {
  if (curve != curvecopy) {
    mont_pt_copy(p, &curve->P, &curvecopy->P);
    mont_pt_copy(p, &curve->Q, &curvecopy->Q);
    fp2_Copy(p, &curve->a, &curvecopy->a);
    fp2_Copy(p, &curve->b, &curvecopy->b);
  }
}

void mont_curve_init(ff_Params* p, mont_curve_int_t* curve) {
  curve->ffData = p;

  mont_pt_init(p, &curve->P);
  mont_pt_init(p, &curve->Q);
  fp2_Init(p, &curve->a);
  fp2_Init(p, &curve->b);
}

void mp_export(void *rop, size_t *countp, int order, size_t size, int endian, size_t nails, const mpz_t op) {
  mpz_export(rop, countp, order, size, endian, nails, op);
}

static void itoos(const mp to_enc, unsigned char* enc) {
  mp_export(enc, NULL, -1, 1, 1, 0, to_enc);
}


static void fptoos(const mp to_enc, unsigned char* enc) {
  itoos(to_enc, enc);
}

static size_t get_np_len(const mp p) {
  return BITS_TO_BYTES_CEIL(mp_sizeinbase(p, 2));
}


int main()
{
	const sike_params_raw_t *params_raw = NULL;
	params_raw = &SIKEp751;
	sike_params_t params;
	sike_setup_params(params_raw, &params);
  mont_curve_int_t* EA = &params.EA;
	//int rc =0;
	//rc = test_sike_int(params_raw->name, &params);
/*
printf("To check if things work properly \n");
printf("The points P and Q of of type mont_pt_t in EA \n");
printf("P (mont_pt_t)= \n");
  gmp_printf(" %Zx \n",(EA->P).x.x0);
  printf("\n");
  printf("Q (mont_pt_t)= \n");
  gmp_printf(" %Zx \n",(EA->Q).x.x0);
printf("\n");
  printf("a (fp2)= \n");
  gmp_printf(" %Zx \n",(EA->a).x1);
printf("\n");
  printf("b (fp2)= \n");
  gmp_printf(" %Zx \n",(EA->b).x1);
*/
  printf("initalize x0");
  mont_pt_t T;
  mont_pt_t P,Q;
  mpz_init(T.x.x0);
  mpz_init(T.x.x1);
  mpz_init(T.y.x0);
  mpz_init(T.y.x1);
  mpz_init(P.x.x0);
	mpz_init(P.x.x1);
	mpz_init(P.y.x0);
	mpz_init(P.y.x1);
  mpz_init(Q.x.x0);
	mpz_init(Q.x.x1);
	mpz_init(Q.y.x0);
	mpz_init(Q.y.x1);
  printf("Value Initalized");
  mpz_set_str((T.x.x0),"0x77B3BB69009428A327D43CA60169715F547454F88CD213452DF58A7252C2B3C3D00D52CCD3133D54041D8BCAEA291F2057202328712CD395575CD7CCD3CE70C0A1EBF633BA946559458878F41F9FDD1727E2C31125B2FE5B713067093331",0);
  mpz_set_str((T.x.x1),"0x6D91393A57DBF47FD6DCF841F17ECD719CAE1D33C6832A75B0F168855BCC38D2A4792DFF9BC86DEACA10B1AA808D539B167D73BBA32168687FA3F85AE93A1ADDE5BD1FD5B681DCC6C34454D4496976C22D80C95E42B12576FC0FB4074B9F",0);
	mpz_set_str((T.y.x0),"0x53B55053E3F04FC315EFB1B7B2C4AFCB4FEF12CE744AF3B243C6E6B1417E94A78D4980DDE181896464923E01AACC3DA040A0747CA67554A352684DA207C49022D930732DF6BD0BF37E1F5C16917669A70F88059C1C739A79D7CFA0C529D9",0);
	mpz_set_str((T.y.x1),"0x44E44196909252ECD7B9164323815294F02AED22C4E4EB43D2CE2BC5F29EB575D45CA8B6B4C4242E369AE3A1EFC844E9D1C57B0AE3374BC2CEDAD16B0C699158332E2D9AB3F0025C0348C5F70FDC4DD7C4865E64B8B843F03D807447D5EF",0);

  mpz_set_str((P.x.x0),"0x77B3BB69009428A327D43CA60169715F547454F88CD213452DF58A7252C2B3C3D00D52CCD3133D54041D8BCAEA291F2057202328712CD395575CD7CCD3CE70C0A1EBF633BA946559458878F41F9FDD1727E2C31125B2FE5B713067093331",0);
  mpz_set_str((P.x.x1),"0x6D91393A57DBF47FD6DCF841F17ECD719CAE1D33C6832A75B0F168855BCC38D2A4792DFF9BC86DEACA10B1AA808D539B167D73BBA32168687FA3F85AE93A1ADDE5BD1FD5B681DCC6C34454D4496976C22D80C95E42B12576FC0FB4074B9F",0);
	mpz_set_str((P.y.x0),"0x53B55053E3F04FC315EFB1B7B2C4AFCB4FEF12CE744AF3B243C6E6B1417E94A78D4980DDE181896464923E01AACC3DA040A0747CA67554A352684DA207C49022D930732DF6BD0BF37E1F5C16917669A70F88059C1C739A79D7CFA0C529D9",0);
	mpz_set_str((P.y.x1),"0x44E44196909252ECD7B9164323815294F02AED22C4E4EB43D2CE2BC5F29EB575D45CA8B6B4C4242E369AE3A1EFC844E9D1C57B0AE3374BC2CEDAD16B0C699158332E2D9AB3F0025C0348C5F70FDC4DD7C4865E64B8B843F03D807447D5EF",0);


  for(int i=0;i<100;i++) {
    xDBL(EA,&T,&T);
    gmp_printf("T > %Zx \n",T);
  }

  for(int i=0;i<100;i++) {
    xADD(EA,&P,&T,&Q);
  }

  for(int i=0;i<100;i++) {
    xTPL(EA,&P,&Q);
  }


	return 0;
}




