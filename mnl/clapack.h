#ifndef __CLAPACK_H
#define __CLAPACK_H

extern "C" {
#include "f2c.h"
 
/* Subroutine */ int LAPACKMANGLE(cbdsqr)(char *uplo, integer *n, integer *ncvt, integer *
	nru, integer *ncc, real *d__, real *e, complex *vt, integer *ldvt, 
	complex *u, integer *ldu, complex *c__, integer *ldc, real *rwork, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cgbbrd)(char *vect, integer *m, integer *n, integer *ncc,
	 integer *kl, integer *ku, complex *ab, integer *ldab, real *d__, 
	real *e, complex *q, integer *ldq, complex *pt, integer *ldpt, 
	complex *c__, integer *ldc, complex *work, real *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cgbcon)(char *norm, integer *n, integer *kl, integer *ku,
	 complex *ab, integer *ldab, integer *ipiv, real *anorm, real *rcond, 
	complex *work, real *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cgbequ)(integer *m, integer *n, integer *kl, integer *ku,
	 complex *ab, integer *ldab, real *r__, real *c__, real *rowcnd, real 
	*colcnd, real *amax, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cgbrfs)(char *trans, integer *n, integer *kl, integer *
	ku, integer *nrhs, complex *ab, integer *ldab, complex *afb, integer *
	ldafb, integer *ipiv, complex *b, integer *ldb, complex *x, integer *
	ldx, real *ferr, real *berr, complex *work, real *rwork, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(cgbsv)(integer *n, integer *kl, integer *ku, integer *
	nrhs, complex *ab, integer *ldab, integer *ipiv, complex *b, integer *
	ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cgbsvx)(char *fact, char *trans, integer *n, integer *kl,
	 integer *ku, integer *nrhs, complex *ab, integer *ldab, complex *afb,
	 integer *ldafb, integer *ipiv, char *equed, real *r__, real *c__, 
	complex *b, integer *ldb, complex *x, integer *ldx, real *rcond, real 
	*ferr, real *berr, complex *work, real *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cgbtf2)(integer *m, integer *n, integer *kl, integer *ku,
	 complex *ab, integer *ldab, integer *ipiv, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cgbtrf)(integer *m, integer *n, integer *kl, integer *ku,
	 complex *ab, integer *ldab, integer *ipiv, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cgbtrs)(char *trans, integer *n, integer *kl, integer *
	ku, integer *nrhs, complex *ab, integer *ldab, integer *ipiv, complex 
	*b, integer *ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cgebak)(char *job, char *side, integer *n, integer *ilo, 
	integer *ihi, real *scale, integer *m, complex *v, integer *ldv, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cgebal)(char *job, integer *n, complex *a, integer *lda, 
	integer *ilo, integer *ihi, real *scale, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cgebd2)(integer *m, integer *n, complex *a, integer *lda,
	 real *d__, real *e, complex *tauq, complex *taup, complex *work, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cgebrd)(integer *m, integer *n, complex *a, integer *lda,
	 real *d__, real *e, complex *tauq, complex *taup, complex *work, 
	integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cgecon)(char *norm, integer *n, complex *a, integer *lda,
	 real *anorm, real *rcond, complex *work, real *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cgeequ)(integer *m, integer *n, complex *a, integer *lda,
	 real *r__, real *c__, real *rowcnd, real *colcnd, real *amax, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cgees)(char *jobvs, char *sort, L_fp select, integer *n, 
	complex *a, integer *lda, integer *sdim, complex *w, complex *vs, 
	integer *ldvs, complex *work, integer *lwork, real *rwork, logical *
	bwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cgeesx)(char *jobvs, char *sort, L_fp select, char *
	sense, integer *n, complex *a, integer *lda, integer *sdim, complex *
	w, complex *vs, integer *ldvs, real *rconde, real *rcondv, complex *
	work, integer *lwork, real *rwork, logical *bwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cgeev)(char *jobvl, char *jobvr, integer *n, complex *a, 
	integer *lda, complex *w, complex *vl, integer *ldvl, complex *vr, 
	integer *ldvr, complex *work, integer *lwork, real *rwork, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(cgeevx)(char *balanc, char *jobvl, char *jobvr, char *
	sense, integer *n, complex *a, integer *lda, complex *w, complex *vl, 
	integer *ldvl, complex *vr, integer *ldvr, integer *ilo, integer *ihi,
	 real *scale, real *abnrm, real *rconde, real *rcondv, complex *work, 
	integer *lwork, real *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cgegs)(char *jobvsl, char *jobvsr, integer *n, complex *
	a, integer *lda, complex *b, integer *ldb, complex *alpha, complex *
	beta, complex *vsl, integer *ldvsl, complex *vsr, integer *ldvsr, 
	complex *work, integer *lwork, real *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cgegv)(char *jobvl, char *jobvr, integer *n, complex *a, 
	integer *lda, complex *b, integer *ldb, complex *alpha, complex *beta,
	 complex *vl, integer *ldvl, complex *vr, integer *ldvr, complex *
	work, integer *lwork, real *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cgehd2)(integer *n, integer *ilo, integer *ihi, complex *
	a, integer *lda, complex *tau, complex *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cgehrd)(integer *n, integer *ilo, integer *ihi, complex *
	a, integer *lda, complex *tau, complex *work, integer *lwork, integer 
	*info);
 
/* Subroutine */ int LAPACKMANGLE(cgelq2)(integer *m, integer *n, complex *a, integer *lda,
	 complex *tau, complex *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cgelqf)(integer *m, integer *n, complex *a, integer *lda,
	 complex *tau, complex *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cgels)(char *trans, integer *m, integer *n, integer *
	nrhs, complex *a, integer *lda, complex *b, integer *ldb, complex *
	work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cgelsx)(integer *m, integer *n, integer *nrhs, complex *
	a, integer *lda, complex *b, integer *ldb, integer *jpvt, real *rcond,
	 integer *rank, complex *work, real *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cgelsy)(integer *m, integer *n, integer *nrhs, complex *
	a, integer *lda, complex *b, integer *ldb, integer *jpvt, real *rcond,
	 integer *rank, complex *work, integer *lwork, real *rwork, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(cgeql2)(integer *m, integer *n, complex *a, integer *lda,
	 complex *tau, complex *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cgeqlf)(integer *m, integer *n, complex *a, integer *lda,
	 complex *tau, complex *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cgeqp3)(integer *m, integer *n, complex *a, integer *lda,
	 integer *jpvt, complex *tau, complex *work, integer *lwork, real *
	rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cgeqpf)(integer *m, integer *n, complex *a, integer *lda,
	 integer *jpvt, complex *tau, complex *work, real *rwork, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(cgeqr2)(integer *m, integer *n, complex *a, integer *lda,
	 complex *tau, complex *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cgeqrf)(integer *m, integer *n, complex *a, integer *lda,
	 complex *tau, complex *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cgerfs)(char *trans, integer *n, integer *nrhs, complex *
	a, integer *lda, complex *af, integer *ldaf, integer *ipiv, complex *
	b, integer *ldb, complex *x, integer *ldx, real *ferr, real *berr, 
	complex *work, real *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cgerq2)(integer *m, integer *n, complex *a, integer *lda,
	 complex *tau, complex *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cgerqf)(integer *m, integer *n, complex *a, integer *lda,
	 complex *tau, complex *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cgesc2)(integer *n, complex *a, integer *lda, complex *
	rhs, integer *ipiv, integer *jpiv, real *scale);
 
/* Subroutine */ int LAPACKMANGLE(cgesv)(integer *n, integer *nrhs, complex *a, integer *
	lda, integer *ipiv, complex *b, integer *ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cgesvx)(char *fact, char *trans, integer *n, integer *
	nrhs, complex *a, integer *lda, complex *af, integer *ldaf, integer *
	ipiv, char *equed, real *r__, real *c__, complex *b, integer *ldb, 
	complex *x, integer *ldx, real *rcond, real *ferr, real *berr, 
	complex *work, real *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cgetc2)(integer *n, complex *a, integer *lda, integer *
	ipiv, integer *jpiv, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cgetf2)(integer *m, integer *n, complex *a, integer *lda,
	 integer *ipiv, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cgetrf)(integer *m, integer *n, complex *a, integer *lda,
	 integer *ipiv, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cgetri)(integer *n, complex *a, integer *lda, integer *
	ipiv, complex *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cgetrs)(char *trans, integer *n, integer *nrhs, complex *
	a, integer *lda, integer *ipiv, complex *b, integer *ldb, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(cggbak)(char *job, char *side, integer *n, integer *ilo, 
	integer *ihi, real *lscale, real *rscale, integer *m, complex *v, 
	integer *ldv, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cggbal)(char *job, integer *n, complex *a, integer *lda, 
	complex *b, integer *ldb, integer *ilo, integer *ihi, real *lscale, 
	real *rscale, real *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cgges)(char *jobvsl, char *jobvsr, char *sort, L_fp 
	selctg, integer *n, complex *a, integer *lda, complex *b, integer *
	ldb, integer *sdim, complex *alpha, complex *beta, complex *vsl, 
	integer *ldvsl, complex *vsr, integer *ldvsr, complex *work, integer *
	lwork, real *rwork, logical *bwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cggesx)(char *jobvsl, char *jobvsr, char *sort, L_fp 
	selctg, char *sense, integer *n, complex *a, integer *lda, complex *b,
	 integer *ldb, integer *sdim, complex *alpha, complex *beta, complex *
	vsl, integer *ldvsl, complex *vsr, integer *ldvsr, real *rconde, real 
	*rcondv, complex *work, integer *lwork, real *rwork, integer *iwork, 
	integer *liwork, logical *bwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cggev)(char *jobvl, char *jobvr, integer *n, complex *a, 
	integer *lda, complex *b, integer *ldb, complex *alpha, complex *beta,
	 complex *vl, integer *ldvl, complex *vr, integer *ldvr, complex *
	work, integer *lwork, real *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cggevx)(char *balanc, char *jobvl, char *jobvr, char *
	sense, integer *n, complex *a, integer *lda, complex *b, integer *ldb,
	 complex *alpha, complex *beta, complex *vl, integer *ldvl, complex *
	vr, integer *ldvr, integer *ilo, integer *ihi, real *lscale, real *
	rscale, real *abnrm, real *bbnrm, real *rconde, real *rcondv, complex 
	*work, integer *lwork, real *rwork, integer *iwork, logical *bwork, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cggglm)(integer *n, integer *m, integer *p, complex *a, 
	integer *lda, complex *b, integer *ldb, complex *d__, complex *x, 
	complex *y, complex *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cgghrd)(char *compq, char *compz, integer *n, integer *
	ilo, integer *ihi, complex *a, integer *lda, complex *b, integer *ldb,
	 complex *q, integer *ldq, complex *z__, integer *ldz, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cgglse)(integer *m, integer *n, integer *p, complex *a, 
	integer *lda, complex *b, integer *ldb, complex *c__, complex *d__, 
	complex *x, complex *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cggqrf)(integer *n, integer *m, integer *p, complex *a, 
	integer *lda, complex *taua, complex *b, integer *ldb, complex *taub, 
	complex *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cggrqf)(integer *m, integer *p, integer *n, complex *a, 
	integer *lda, complex *taua, complex *b, integer *ldb, complex *taub, 
	complex *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cggsvd)(char *jobu, char *jobv, char *jobq, integer *m, 
	integer *n, integer *p, integer *k, integer *l, complex *a, integer *
	lda, complex *b, integer *ldb, real *alpha, real *beta, complex *u, 
	integer *ldu, complex *v, integer *ldv, complex *q, integer *ldq, 
	complex *work, real *rwork, integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cggsvp)(char *jobu, char *jobv, char *jobq, integer *m, 
	integer *p, integer *n, complex *a, integer *lda, complex *b, integer 
	*ldb, real *tola, real *tolb, integer *k, integer *l, complex *u, 
	integer *ldu, complex *v, integer *ldv, complex *q, integer *ldq, 
	integer *iwork, real *rwork, complex *tau, complex *work, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(cgtcon)(char *norm, integer *n, complex *dl, complex *
	d__, complex *du, complex *du2, integer *ipiv, real *anorm, real *
	rcond, complex *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cgtrfs)(char *trans, integer *n, integer *nrhs, complex *
	dl, complex *d__, complex *du, complex *dlf, complex *df, complex *
	duf, complex *du2, integer *ipiv, complex *b, integer *ldb, complex *
	x, integer *ldx, real *ferr, real *berr, complex *work, real *rwork, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cgtsv)(integer *n, integer *nrhs, complex *dl, complex *
	d__, complex *du, complex *b, integer *ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cgtsvx)(char *fact, char *trans, integer *n, integer *
	nrhs, complex *dl, complex *d__, complex *du, complex *dlf, complex *
	df, complex *duf, complex *du2, integer *ipiv, complex *b, integer *
	ldb, complex *x, integer *ldx, real *rcond, real *ferr, real *berr, 
	complex *work, real *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cgttrf)(integer *n, complex *dl, complex *d__, complex *
	du, complex *du2, integer *ipiv, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cgttrs)(char *trans, integer *n, integer *nrhs, complex *
	dl, complex *d__, complex *du, complex *du2, integer *ipiv, complex *
	b, integer *ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cgtts2)(integer *itrans, integer *n, integer *nrhs, 
	complex *dl, complex *d__, complex *du, complex *du2, integer *ipiv, 
	complex *b, integer *ldb);
 
/* Subroutine */ int LAPACKMANGLE(chbev)(char *jobz, char *uplo, integer *n, integer *kd, 
	complex *ab, integer *ldab, real *w, complex *z__, integer *ldz, 
	complex *work, real *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(chbevd)(char *jobz, char *uplo, integer *n, integer *kd, 
	complex *ab, integer *ldab, real *w, complex *z__, integer *ldz, 
	complex *work, integer *lwork, real *rwork, integer *lrwork, integer *
	iwork, integer *liwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(chbevx)(char *jobz, char *range, char *uplo, integer *n, 
	integer *kd, complex *ab, integer *ldab, complex *q, integer *ldq, 
	real *vl, real *vu, integer *il, integer *iu, real *abstol, integer *
	m, real *w, complex *z__, integer *ldz, complex *work, real *rwork, 
	integer *iwork, integer *ifail, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(chbgst)(char *vect, char *uplo, integer *n, integer *ka, 
	integer *kb, complex *ab, integer *ldab, complex *bb, integer *ldbb, 
	complex *x, integer *ldx, complex *work, real *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(chbgv)(char *jobz, char *uplo, integer *n, integer *ka, 
	integer *kb, complex *ab, integer *ldab, complex *bb, integer *ldbb, 
	real *w, complex *z__, integer *ldz, complex *work, real *rwork, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(chbgvx)(char *jobz, char *range, char *uplo, integer *n, 
	integer *ka, integer *kb, complex *ab, integer *ldab, complex *bb, 
	integer *ldbb, complex *q, integer *ldq, real *vl, real *vu, integer *
	il, integer *iu, real *abstol, integer *m, real *w, complex *z__, 
	integer *ldz, complex *work, real *rwork, integer *iwork, integer *
	ifail, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(chbtrd)(char *vect, char *uplo, integer *n, integer *kd, 
	complex *ab, integer *ldab, real *d__, real *e, complex *q, integer *
	ldq, complex *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(checon)(char *uplo, integer *n, complex *a, integer *lda,
	 integer *ipiv, real *anorm, real *rcond, complex *work, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(cheev)(char *jobz, char *uplo, integer *n, complex *a, 
	integer *lda, real *w, complex *work, integer *lwork, real *rwork, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cheevd)(char *jobz, char *uplo, integer *n, complex *a, 
	integer *lda, real *w, complex *work, integer *lwork, real *rwork, 
	integer *lrwork, integer *iwork, integer *liwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cheevr)(char *jobz, char *range, char *uplo, integer *n, 
	complex *a, integer *lda, real *vl, real *vu, integer *il, integer *
	iu, real *abstol, integer *m, real *w, complex *z__, integer *ldz, 
	integer *isuppz, complex *work, integer *lwork, real *rwork, integer *
	lrwork, integer *iwork, integer *liwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cheevx)(char *jobz, char *range, char *uplo, integer *n, 
	complex *a, integer *lda, real *vl, real *vu, integer *il, integer *
	iu, real *abstol, integer *m, real *w, complex *z__, integer *ldz, 
	complex *work, integer *lwork, real *rwork, integer *iwork, integer *
	ifail, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(chegs2)(integer *itype, char *uplo, integer *n, complex *
	a, integer *lda, complex *b, integer *ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(chegst)(integer *itype, char *uplo, integer *n, complex *
	a, integer *lda, complex *b, integer *ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(chegv)(integer *itype, char *jobz, char *uplo, integer *
	n, complex *a, integer *lda, complex *b, integer *ldb, real *w, 
	complex *work, integer *lwork, real *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(chegvd)(integer *itype, char *jobz, char *uplo, integer *
	n, complex *a, integer *lda, complex *b, integer *ldb, real *w, 
	complex *work, integer *lwork, real *rwork, integer *lrwork, integer *
	iwork, integer *liwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(chegvx)(integer *itype, char *jobz, char *range, char *
	uplo, integer *n, complex *a, integer *lda, complex *b, integer *ldb, 
	real *vl, real *vu, integer *il, integer *iu, real *abstol, integer *
	m, real *w, complex *z__, integer *ldz, complex *work, integer *lwork,
	 real *rwork, integer *iwork, integer *ifail, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cherfs)(char *uplo, integer *n, integer *nrhs, complex *
	a, integer *lda, complex *af, integer *ldaf, integer *ipiv, complex *
	b, integer *ldb, complex *x, integer *ldx, real *ferr, real *berr, 
	complex *work, real *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(chesv)(char *uplo, integer *n, integer *nrhs, complex *a,
	 integer *lda, integer *ipiv, complex *b, integer *ldb, complex *work,
	 integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(chesvx)(char *fact, char *uplo, integer *n, integer *
	nrhs, complex *a, integer *lda, complex *af, integer *ldaf, integer *
	ipiv, complex *b, integer *ldb, complex *x, integer *ldx, real *rcond,
	 real *ferr, real *berr, complex *work, integer *lwork, real *rwork, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(chetf2)(char *uplo, integer *n, complex *a, integer *lda,
	 integer *ipiv, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(chetrd)(char *uplo, integer *n, complex *a, integer *lda,
	 real *d__, real *e, complex *tau, complex *work, integer *lwork, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(chetrf)(char *uplo, integer *n, complex *a, integer *lda,
	 integer *ipiv, complex *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(chetri)(char *uplo, integer *n, complex *a, integer *lda,
	 integer *ipiv, complex *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(chetrs)(char *uplo, integer *n, integer *nrhs, complex *
	a, integer *lda, integer *ipiv, complex *b, integer *ldb, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(chgeqz)(char *job, char *compq, char *compz, integer *n, 
	integer *ilo, integer *ihi, complex *a, integer *lda, complex *b, 
	integer *ldb, complex *alpha, complex *beta, complex *q, integer *ldq,
	 complex *z__, integer *ldz, complex *work, integer *lwork, real *
	rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(chpcon)(char *uplo, integer *n, complex *ap, integer *
	ipiv, real *anorm, real *rcond, complex *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(chpev)(char *jobz, char *uplo, integer *n, complex *ap, 
	real *w, complex *z__, integer *ldz, complex *work, real *rwork, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(chpevd)(char *jobz, char *uplo, integer *n, complex *ap, 
	real *w, complex *z__, integer *ldz, complex *work, integer *lwork, 
	real *rwork, integer *lrwork, integer *iwork, integer *liwork, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(chpevx)(char *jobz, char *range, char *uplo, integer *n, 
	complex *ap, real *vl, real *vu, integer *il, integer *iu, real *
	abstol, integer *m, real *w, complex *z__, integer *ldz, complex *
	work, real *rwork, integer *iwork, integer *ifail, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(chpgst)(integer *itype, char *uplo, integer *n, complex *
	ap, complex *bp, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(chpgv)(integer *itype, char *jobz, char *uplo, integer *
	n, complex *ap, complex *bp, real *w, complex *z__, integer *ldz, 
	complex *work, real *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(chpgvd)(integer *itype, char *jobz, char *uplo, integer *
	n, complex *ap, complex *bp, real *w, complex *z__, integer *ldz, 
	complex *work, integer *lwork, real *rwork, integer *lrwork, integer *
	iwork, integer *liwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(chpgvx)(integer *itype, char *jobz, char *range, char *
	uplo, integer *n, complex *ap, complex *bp, real *vl, real *vu, 
	integer *il, integer *iu, real *abstol, integer *m, real *w, complex *
	z__, integer *ldz, complex *work, real *rwork, integer *iwork, 
	integer *ifail, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(chprfs)(char *uplo, integer *n, integer *nrhs, complex *
	ap, complex *afp, integer *ipiv, complex *b, integer *ldb, complex *x,
	 integer *ldx, real *ferr, real *berr, complex *work, real *rwork, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(chpsv)(char *uplo, integer *n, integer *nrhs, complex *
	ap, integer *ipiv, complex *b, integer *ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(chpsvx)(char *fact, char *uplo, integer *n, integer *
	nrhs, complex *ap, complex *afp, integer *ipiv, complex *b, integer *
	ldb, complex *x, integer *ldx, real *rcond, real *ferr, real *berr, 
	complex *work, real *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(chptrd)(char *uplo, integer *n, complex *ap, real *d__, 
	real *e, complex *tau, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(chptrf)(char *uplo, integer *n, complex *ap, integer *
	ipiv, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(chptri)(char *uplo, integer *n, complex *ap, integer *
	ipiv, complex *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(chptrs)(char *uplo, integer *n, integer *nrhs, complex *
	ap, integer *ipiv, complex *b, integer *ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(chsein)(char *side, char *eigsrc, char *initv, logical *
	select, integer *n, complex *h__, integer *ldh, complex *w, complex *
	vl, integer *ldvl, complex *vr, integer *ldvr, integer *mm, integer *
	m, complex *work, real *rwork, integer *ifaill, integer *ifailr, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(chseqr)(char *job, char *compz, integer *n, integer *ilo,
	 integer *ihi, complex *h__, integer *ldh, complex *w, complex *z__, 
	integer *ldz, complex *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(clabrd)(integer *m, integer *n, integer *nb, complex *a, 
	integer *lda, real *d__, real *e, complex *tauq, complex *taup, 
	complex *x, integer *ldx, complex *y, integer *ldy);
 
/* Subroutine */ int LAPACKMANGLE(clacgv)(integer *n, complex *x, integer *incx);
 
/* Subroutine */ int LAPACKMANGLE(clacon)(integer *n, complex *v, complex *x, real *est, 
	integer *kase);
 
/* Subroutine */ int LAPACKMANGLE(clacp2)(char *uplo, integer *m, integer *n, real *a, 
	integer *lda, complex *b, integer *ldb);
 
/* Subroutine */ int LAPACKMANGLE(clacpy)(char *uplo, integer *m, integer *n, complex *a, 
	integer *lda, complex *b, integer *ldb);
 
/* Subroutine */ int LAPACKMANGLE(clacrm)(integer *m, integer *n, complex *a, integer *lda,
	 real *b, integer *ldb, complex *c__, integer *ldc, real *rwork);
 
/* Subroutine */ int LAPACKMANGLE(clacrt)(integer *n, complex *cx, integer *incx, complex *
	cy, integer *incy, complex *c__, complex *s);
 
/* Subroutine */ int LAPACKMANGLE(claed0)(integer *qsiz, integer *n, real *d__, real *e, 
	complex *q, integer *ldq, complex *qstore, integer *ldqs, real *rwork,
	 integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(claed7)(integer *n, integer *cutpnt, integer *qsiz, 
	integer *tlvls, integer *curlvl, integer *curpbm, real *d__, complex *
	q, integer *ldq, real *rho, integer *indxq, real *qstore, integer *
	qptr, integer *prmptr, integer *perm, integer *givptr, integer *
	givcol, real *givnum, complex *work, real *rwork, integer *iwork, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(claed8)(integer *k, integer *n, integer *qsiz, complex *
	q, integer *ldq, real *d__, real *rho, integer *cutpnt, real *z__, 
	real *dlamda, complex *q2, integer *ldq2, real *w, integer *indxp, 
	integer *indx, integer *indxq, integer *perm, integer *givptr, 
	integer *givcol, real *givnum, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(claein)(logical *rightv, logical *noinit, integer *n, 
	complex *h__, integer *ldh, complex *w, complex *v, complex *b, 
	integer *ldb, real *rwork, real *eps3, real *smlnum, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(claesy)(complex *a, complex *b, complex *c__, complex *
	rt1, complex *rt2, complex *evscal, complex *cs1, complex *sn1);
 
/* Subroutine */ int LAPACKMANGLE(claev2)(complex *a, complex *b, complex *c__, real *rt1, 
	real *rt2, real *cs1, complex *sn1);
 
/* Subroutine */ int LAPACKMANGLE(clags2)(logical *upper, real *a1, complex *a2, real *a3, 
	real *b1, complex *b2, real *b3, real *csu, complex *snu, real *csv, 
	complex *snv, real *csq, complex *snq);
 
/* Subroutine */ int LAPACKMANGLE(clagtm)(char *trans, integer *n, integer *nrhs, real *
	alpha, complex *dl, complex *d__, complex *du, complex *x, integer *
	ldx, real *beta, complex *b, integer *ldb);
 
/* Subroutine */ int LAPACKMANGLE(clahef)(char *uplo, integer *n, integer *nb, integer *kb,
	 complex *a, integer *lda, integer *ipiv, complex *w, integer *ldw, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(clahqr)(logical *wantt, logical *wantz, integer *n, 
	integer *ilo, integer *ihi, complex *h__, integer *ldh, complex *w, 
	integer *iloz, integer *ihiz, complex *z__, integer *ldz, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(clahrd)(integer *n, integer *k, integer *nb, complex *a, 
	integer *lda, complex *tau, complex *t, integer *ldt, complex *y, 
	integer *ldy);
 
/* Subroutine */ int LAPACKMANGLE(claic1)(integer *job, integer *j, complex *x, real *sest,
	 complex *w, complex *gamma, real *sestpr, complex *s, complex *c__);
 
/* Subroutine */ int LAPACKMANGLE(clals0)(integer *icompq, integer *nl, integer *nr, 
	integer *sqre, integer *nrhs, complex *b, integer *ldb, complex *bx, 
	integer *ldbx, integer *perm, integer *givptr, integer *givcol, 
	integer *ldgcol, real *givnum, integer *ldgnum, real *poles, real *
	difl, real *difr, real *z__, integer *k, real *c__, real *s, real *
	rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(clalsa)(integer *icompq, integer *smlsiz, integer *n, 
	integer *nrhs, complex *b, integer *ldb, complex *bx, integer *ldbx, 
	real *u, integer *ldu, real *vt, integer *k, real *difl, real *difr, 
	real *z__, real *poles, integer *givptr, integer *givcol, integer *
	ldgcol, integer *perm, real *givnum, real *c__, real *s, real *rwork, 
	integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(clapll)(integer *n, complex *x, integer *incx, complex *
	y, integer *incy, real *ssmin);
 
/* Subroutine */ int LAPACKMANGLE(clapmt)(logical *forwrd, integer *m, integer *n, complex 
	*x, integer *ldx, integer *k);
 
/* Subroutine */ int LAPACKMANGLE(claqgb)(integer *m, integer *n, integer *kl, integer *ku,
	 complex *ab, integer *ldab, real *r__, real *c__, real *rowcnd, real 
	*colcnd, real *amax, char *equed);
 
/* Subroutine */ int LAPACKMANGLE(claqge)(integer *m, integer *n, complex *a, integer *lda,
	 real *r__, real *c__, real *rowcnd, real *colcnd, real *amax, char *
	equed);
 
/* Subroutine */ int LAPACKMANGLE(claqhb)(char *uplo, integer *n, integer *kd, complex *ab,
	 integer *ldab, real *s, real *scond, real *amax, char *equed);
 
/* Subroutine */ int LAPACKMANGLE(claqhe)(char *uplo, integer *n, complex *a, integer *lda,
	 real *s, real *scond, real *amax, char *equed);
 
/* Subroutine */ int LAPACKMANGLE(claqhp)(char *uplo, integer *n, complex *ap, real *s, 
	real *scond, real *amax, char *equed);
 
/* Subroutine */ int LAPACKMANGLE(claqp2)(integer *m, integer *n, integer *offset, complex 
	*a, integer *lda, integer *jpvt, complex *tau, real *vn1, real *vn2, 
	complex *work);
 
/* Subroutine */ int LAPACKMANGLE(claqps)(integer *m, integer *n, integer *offset, integer 
	*nb, integer *kb, complex *a, integer *lda, integer *jpvt, complex *
	tau, real *vn1, real *vn2, complex *auxv, complex *f, integer *ldf);
 
/* Subroutine */ int LAPACKMANGLE(claqsb)(char *uplo, integer *n, integer *kd, complex *ab,
	 integer *ldab, real *s, real *scond, real *amax, char *equed);
 
/* Subroutine */ int LAPACKMANGLE(claqsp)(char *uplo, integer *n, complex *ap, real *s, 
	real *scond, real *amax, char *equed);
 
/* Subroutine */ int LAPACKMANGLE(claqsy)(char *uplo, integer *n, complex *a, integer *lda,
	 real *s, real *scond, real *amax, char *equed);
 
/* Subroutine */ int LAPACKMANGLE(clar1v)(integer *n, integer *b1, integer *bn, real *
	sigma, real *d__, real *l, real *ld, real *lld, real *gersch, complex 
	*z__, real *ztz, real *mingma, integer *r__, integer *isuppz, real *
	work);
 
/* Subroutine */ int LAPACKMANGLE(clar2v)(integer *n, complex *x, complex *y, complex *z__,
	 integer *incx, real *c__, complex *s, integer *incc);
 
/* Subroutine */ int LAPACKMANGLE(clarcm)(integer *m, integer *n, real *a, integer *lda, 
	complex *b, integer *ldb, complex *c__, integer *ldc, real *rwork);
 
/* Subroutine */ int LAPACKMANGLE(clarf)(char *side, integer *m, integer *n, complex *v, 
	integer *incv, complex *tau, complex *c__, integer *ldc, complex *
	work);
 
/* Subroutine */ int LAPACKMANGLE(clarfb)(char *side, char *trans, char *direct, char *
	storev, integer *m, integer *n, integer *k, complex *v, integer *ldv, 
	complex *t, integer *ldt, complex *c__, integer *ldc, complex *work, 
	integer *ldwork);
 
/* Subroutine */ int LAPACKMANGLE(clarfg)(integer *n, complex *alpha, complex *x, integer *
	incx, complex *tau);
 
/* Subroutine */ int LAPACKMANGLE(clarft)(char *direct, char *storev, integer *n, integer *
	k, complex *v, integer *ldv, complex *tau, complex *t, integer *ldt);
 
/* Subroutine */ int LAPACKMANGLE(clarfx)(char *side, integer *m, integer *n, complex *v, 
	complex *tau, complex *c__, integer *ldc, complex *work);
 
/* Subroutine */ int LAPACKMANGLE(clargv)(integer *n, complex *x, integer *incx, complex *
	y, integer *incy, real *c__, integer *incc);
 
/* Subroutine */ int LAPACKMANGLE(clarnv)(integer *idist, integer *iseed, integer *n, 
	complex *x);
 
/* Subroutine */ int LAPACKMANGLE(clarrv)(integer *n, real *d__, real *l, integer *isplit, 
	integer *m, real *w, integer *iblock, real *gersch, real *tol, 
	complex *z__, integer *ldz, integer *isuppz, real *work, integer *
	iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(clartg)(complex *f, complex *g, real *cs, complex *sn, 
	complex *r__);
 
/* Subroutine */ int LAPACKMANGLE(clartv)(integer *n, complex *x, integer *incx, complex *
	y, integer *incy, real *c__, complex *s, integer *incc);
 
/* Subroutine */ int LAPACKMANGLE(clarz)(char *side, integer *m, integer *n, integer *l, 
	complex *v, integer *incv, complex *tau, complex *c__, integer *ldc, 
	complex *work);
 
/* Subroutine */ int LAPACKMANGLE(clarzb)(char *side, char *trans, char *direct, char *
	storev, integer *m, integer *n, integer *k, integer *l, complex *v, 
	integer *ldv, complex *t, integer *ldt, complex *c__, integer *ldc, 
	complex *work, integer *ldwork);
 
/* Subroutine */ int LAPACKMANGLE(clarzt)(char *direct, char *storev, integer *n, integer *
	k, complex *v, integer *ldv, complex *tau, complex *t, integer *ldt);
 
/* Subroutine */ int LAPACKMANGLE(clascl)(char *type__, integer *kl, integer *ku, real *
	cfrom, real *cto, integer *m, integer *n, complex *a, integer *lda, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(claset)(char *uplo, integer *m, integer *n, complex *
	alpha, complex *beta, complex *a, integer *lda);
 
/* Subroutine */ int LAPACKMANGLE(clasr)(char *side, char *pivot, char *direct, integer *m,
	 integer *n, real *c__, real *s, complex *a, integer *lda);
 
/* Subroutine */ int LAPACKMANGLE(classq)(integer *n, complex *x, integer *incx, real *
	scale, real *sumsq);
 
/* Subroutine */ int LAPACKMANGLE(claswp)(integer *n, complex *a, integer *lda, integer *
	k1, integer *k2, integer *ipiv, integer *incx);
 
/* Subroutine */ int LAPACKMANGLE(clasyf)(char *uplo, integer *n, integer *nb, integer *kb,
	 complex *a, integer *lda, integer *ipiv, complex *w, integer *ldw, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(clatbs)(char *uplo, char *trans, char *diag, char *
	normin, integer *n, integer *kd, complex *ab, integer *ldab, complex *
	x, real *scale, real *cnorm, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(clatdf)(integer *ijob, integer *n, complex *z__, integer 
	*ldz, complex *rhs, real *rdsum, real *rdscal, integer *ipiv, integer 
	*jpiv);
 
/* Subroutine */ int LAPACKMANGLE(clatps)(char *uplo, char *trans, char *diag, char *
	normin, integer *n, complex *ap, complex *x, real *scale, real *cnorm,
	 integer *info);
 
/* Subroutine */ int LAPACKMANGLE(clatrd)(char *uplo, integer *n, integer *nb, complex *a, 
	integer *lda, real *e, complex *tau, complex *w, integer *ldw);
 
/* Subroutine */ int LAPACKMANGLE(clatrs)(char *uplo, char *trans, char *diag, char *
	normin, integer *n, complex *a, integer *lda, complex *x, real *scale,
	 real *cnorm, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(clatrz)(integer *m, integer *n, integer *l, complex *a, 
	integer *lda, complex *tau, complex *work);
 
/* Subroutine */ int LAPACKMANGLE(clatzm)(char *side, integer *m, integer *n, complex *v, 
	integer *incv, complex *tau, complex *c1, complex *c2, integer *ldc, 
	complex *work);
 
/* Subroutine */ int LAPACKMANGLE(clauu2)(char *uplo, integer *n, complex *a, integer *lda,
	 integer *info);
 
/* Subroutine */ int LAPACKMANGLE(clauum)(char *uplo, integer *n, complex *a, integer *lda,
	 integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cpbcon)(char *uplo, integer *n, integer *kd, complex *ab,
	 integer *ldab, real *anorm, real *rcond, complex *work, real *rwork, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cpbequ)(char *uplo, integer *n, integer *kd, complex *ab,
	 integer *ldab, real *s, real *scond, real *amax, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cpbrfs)(char *uplo, integer *n, integer *kd, integer *
	nrhs, complex *ab, integer *ldab, complex *afb, integer *ldafb, 
	complex *b, integer *ldb, complex *x, integer *ldx, real *ferr, real *
	berr, complex *work, real *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cpbstf)(char *uplo, integer *n, integer *kd, complex *ab,
	 integer *ldab, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cpbsv)(char *uplo, integer *n, integer *kd, integer *
	nrhs, complex *ab, integer *ldab, complex *b, integer *ldb, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(cpbsvx)(char *fact, char *uplo, integer *n, integer *kd, 
	integer *nrhs, complex *ab, integer *ldab, complex *afb, integer *
	ldafb, char *equed, real *s, complex *b, integer *ldb, complex *x, 
	integer *ldx, real *rcond, real *ferr, real *berr, complex *work, 
	real *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cpbtf2)(char *uplo, integer *n, integer *kd, complex *ab,
	 integer *ldab, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cpbtrf)(char *uplo, integer *n, integer *kd, complex *ab,
	 integer *ldab, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cpbtrs)(char *uplo, integer *n, integer *kd, integer *
	nrhs, complex *ab, integer *ldab, complex *b, integer *ldb, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(cpocon)(char *uplo, integer *n, complex *a, integer *lda,
	 real *anorm, real *rcond, complex *work, real *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cpoequ)(integer *n, complex *a, integer *lda, real *s, 
	real *scond, real *amax, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cporfs)(char *uplo, integer *n, integer *nrhs, complex *
	a, integer *lda, complex *af, integer *ldaf, complex *b, integer *ldb,
	 complex *x, integer *ldx, real *ferr, real *berr, complex *work, 
	real *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cposv)(char *uplo, integer *n, integer *nrhs, complex *a,
	 integer *lda, complex *b, integer *ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cposvx)(char *fact, char *uplo, integer *n, integer *
	nrhs, complex *a, integer *lda, complex *af, integer *ldaf, char *
	equed, real *s, complex *b, integer *ldb, complex *x, integer *ldx, 
	real *rcond, real *ferr, real *berr, complex *work, real *rwork, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cpotf2)(char *uplo, integer *n, complex *a, integer *lda,
	 integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cpotrf)(char *uplo, integer *n, complex *a, integer *lda,
	 integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cpotri)(char *uplo, integer *n, complex *a, integer *lda,
	 integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cpotrs)(char *uplo, integer *n, integer *nrhs, complex *
	a, integer *lda, complex *b, integer *ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cppcon)(char *uplo, integer *n, complex *ap, real *anorm,
	 real *rcond, complex *work, real *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cppequ)(char *uplo, integer *n, complex *ap, real *s, 
	real *scond, real *amax, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cpprfs)(char *uplo, integer *n, integer *nrhs, complex *
	ap, complex *afp, complex *b, integer *ldb, complex *x, integer *ldx, 
	real *ferr, real *berr, complex *work, real *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cppsv)(char *uplo, integer *n, integer *nrhs, complex *
	ap, complex *b, integer *ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cppsvx)(char *fact, char *uplo, integer *n, integer *
	nrhs, complex *ap, complex *afp, char *equed, real *s, complex *b, 
	integer *ldb, complex *x, integer *ldx, real *rcond, real *ferr, real 
	*berr, complex *work, real *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cpptrf)(char *uplo, integer *n, complex *ap, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(cpptri)(char *uplo, integer *n, complex *ap, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(cpptrs)(char *uplo, integer *n, integer *nrhs, complex *
	ap, complex *b, integer *ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cptcon)(integer *n, real *d__, complex *e, real *anorm, 
	real *rcond, real *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cptrfs)(char *uplo, integer *n, integer *nrhs, real *d__,
	 complex *e, real *df, complex *ef, complex *b, integer *ldb, complex 
	*x, integer *ldx, real *ferr, real *berr, complex *work, real *rwork, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cptsv)(integer *n, integer *nrhs, real *d__, complex *e, 
	complex *b, integer *ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cptsvx)(char *fact, integer *n, integer *nrhs, real *d__,
	 complex *e, real *df, complex *ef, complex *b, integer *ldb, complex 
	*x, integer *ldx, real *rcond, real *ferr, real *berr, complex *work, 
	real *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cpttrf)(integer *n, real *d__, complex *e, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cpttrs)(char *uplo, integer *n, integer *nrhs, real *d__,
	 complex *e, complex *b, integer *ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cptts2)(integer *iuplo, integer *n, integer *nrhs, real *
	d__, complex *e, complex *b, integer *ldb);
 
/* Subroutine */ int LAPACKMANGLE(crot)(integer *n, complex *cx, integer *incx, complex *
	cy, integer *incy, real *c__, complex *s);
 
/* Subroutine */ int LAPACKMANGLE(cspcon)(char *uplo, integer *n, complex *ap, integer *
	ipiv, real *anorm, real *rcond, complex *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cspmv)(char *uplo, integer *n, complex *alpha, complex *
	ap, complex *x, integer *incx, complex *beta, complex *y, integer *
	incy);
 
/* Subroutine */ int LAPACKMANGLE(cspr)(char *uplo, integer *n, complex *alpha, complex *x,
	 integer *incx, complex *ap);
 
/* Subroutine */ int LAPACKMANGLE(csprfs)(char *uplo, integer *n, integer *nrhs, complex *
	ap, complex *afp, integer *ipiv, complex *b, integer *ldb, complex *x,
	 integer *ldx, real *ferr, real *berr, complex *work, real *rwork, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cspsv)(char *uplo, integer *n, integer *nrhs, complex *
	ap, integer *ipiv, complex *b, integer *ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cspsvx)(char *fact, char *uplo, integer *n, integer *
	nrhs, complex *ap, complex *afp, integer *ipiv, complex *b, integer *
	ldb, complex *x, integer *ldx, real *rcond, real *ferr, real *berr, 
	complex *work, real *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(csptrf)(char *uplo, integer *n, complex *ap, integer *
	ipiv, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(csptri)(char *uplo, integer *n, complex *ap, integer *
	ipiv, complex *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(csptrs)(char *uplo, integer *n, integer *nrhs, complex *
	ap, integer *ipiv, complex *b, integer *ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(csrot)(integer *n, complex *cx, integer *incx, complex *
	cy, integer *incy, real *c__, real *s);
 
/* Subroutine */ int LAPACKMANGLE(csrscl)(integer *n, real *sa, complex *sx, integer *incx);
 
/* Subroutine */ int LAPACKMANGLE(cstedc)(char *compz, integer *n, real *d__, real *e, 
	complex *z__, integer *ldz, complex *work, integer *lwork, real *
	rwork, integer *lrwork, integer *iwork, integer *liwork, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(cstein)(integer *n, real *d__, real *e, integer *m, real 
	*w, integer *iblock, integer *isplit, complex *z__, integer *ldz, 
	real *work, integer *iwork, integer *ifail, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(csteqr)(char *compz, integer *n, real *d__, real *e, 
	complex *z__, integer *ldz, real *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(csycon)(char *uplo, integer *n, complex *a, integer *lda,
	 integer *ipiv, real *anorm, real *rcond, complex *work, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(csymv)(char *uplo, integer *n, complex *alpha, complex *
	a, integer *lda, complex *x, integer *incx, complex *beta, complex *y,
	 integer *incy);
 
/* Subroutine */ int LAPACKMANGLE(csyr)(char *uplo, integer *n, complex *alpha, complex *x,
	 integer *incx, complex *a, integer *lda);
 
/* Subroutine */ int LAPACKMANGLE(csyrfs)(char *uplo, integer *n, integer *nrhs, complex *
	a, integer *lda, complex *af, integer *ldaf, integer *ipiv, complex *
	b, integer *ldb, complex *x, integer *ldx, real *ferr, real *berr, 
	complex *work, real *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(csysv)(char *uplo, integer *n, integer *nrhs, complex *a,
	 integer *lda, integer *ipiv, complex *b, integer *ldb, complex *work,
	 integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(csysvx)(char *fact, char *uplo, integer *n, integer *
	nrhs, complex *a, integer *lda, complex *af, integer *ldaf, integer *
	ipiv, complex *b, integer *ldb, complex *x, integer *ldx, real *rcond,
	 real *ferr, real *berr, complex *work, integer *lwork, real *rwork, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(csytf2)(char *uplo, integer *n, complex *a, integer *lda,
	 integer *ipiv, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(csytrf)(char *uplo, integer *n, complex *a, integer *lda,
	 integer *ipiv, complex *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(csytri)(char *uplo, integer *n, complex *a, integer *lda,
	 integer *ipiv, complex *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(csytrs)(char *uplo, integer *n, integer *nrhs, complex *
	a, integer *lda, integer *ipiv, complex *b, integer *ldb, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(ctbcon)(char *norm, char *uplo, char *diag, integer *n, 
	integer *kd, complex *ab, integer *ldab, real *rcond, complex *work, 
	real *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ctbrfs)(char *uplo, char *trans, char *diag, integer *n, 
	integer *kd, integer *nrhs, complex *ab, integer *ldab, complex *b, 
	integer *ldb, complex *x, integer *ldx, real *ferr, real *berr, 
	complex *work, real *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ctbtrs)(char *uplo, char *trans, char *diag, integer *n, 
	integer *kd, integer *nrhs, complex *ab, integer *ldab, complex *b, 
	integer *ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ctgevc)(char *side, char *howmny, logical *select, 
	integer *n, complex *a, integer *lda, complex *b, integer *ldb, 
	complex *vl, integer *ldvl, complex *vr, integer *ldvr, integer *mm, 
	integer *m, complex *work, real *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ctgex2)(logical *wantq, logical *wantz, integer *n, 
	complex *a, integer *lda, complex *b, integer *ldb, complex *q, 
	integer *ldq, complex *z__, integer *ldz, integer *j1, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ctgexc)(logical *wantq, logical *wantz, integer *n, 
	complex *a, integer *lda, complex *b, integer *ldb, complex *q, 
	integer *ldq, complex *z__, integer *ldz, integer *ifst, integer *
	ilst, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ctgsen)(integer *ijob, logical *wantq, logical *wantz, 
	logical *select, integer *n, complex *a, integer *lda, complex *b, 
	integer *ldb, complex *alpha, complex *beta, complex *q, integer *ldq,
	 complex *z__, integer *ldz, integer *m, real *pl, real *pr, real *
	dif, complex *work, integer *lwork, integer *iwork, integer *liwork, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ctgsja)(char *jobu, char *jobv, char *jobq, integer *m, 
	integer *p, integer *n, integer *k, integer *l, complex *a, integer *
	lda, complex *b, integer *ldb, real *tola, real *tolb, real *alpha, 
	real *beta, complex *u, integer *ldu, complex *v, integer *ldv, 
	complex *q, integer *ldq, complex *work, integer *ncycle, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(ctgsna)(char *job, char *howmny, logical *select, 
	integer *n, complex *a, integer *lda, complex *b, integer *ldb, 
	complex *vl, integer *ldvl, complex *vr, integer *ldvr, real *s, real 
	*dif, integer *mm, integer *m, complex *work, integer *lwork, integer 
	*iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ctgsy2)(char *trans, integer *ijob, integer *m, integer *
	n, complex *a, integer *lda, complex *b, integer *ldb, complex *c__, 
	integer *ldc, complex *d__, integer *ldd, complex *e, integer *lde, 
	complex *f, integer *ldf, real *scale, real *rdsum, real *rdscal, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ctgsyl)(char *trans, integer *ijob, integer *m, integer *
	n, complex *a, integer *lda, complex *b, integer *ldb, complex *c__, 
	integer *ldc, complex *d__, integer *ldd, complex *e, integer *lde, 
	complex *f, integer *ldf, real *scale, real *dif, complex *work, 
	integer *lwork, integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ctpcon)(char *norm, char *uplo, char *diag, integer *n, 
	complex *ap, real *rcond, complex *work, real *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ctprfs)(char *uplo, char *trans, char *diag, integer *n, 
	integer *nrhs, complex *ap, complex *b, integer *ldb, complex *x, 
	integer *ldx, real *ferr, real *berr, complex *work, real *rwork, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ctptri)(char *uplo, char *diag, integer *n, complex *ap, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ctptrs)(char *uplo, char *trans, char *diag, integer *n, 
	integer *nrhs, complex *ap, complex *b, integer *ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ctrcon)(char *norm, char *uplo, char *diag, integer *n, 
	complex *a, integer *lda, real *rcond, complex *work, real *rwork, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ctrevc)(char *side, char *howmny, logical *select, 
	integer *n, complex *t, integer *ldt, complex *vl, integer *ldvl, 
	complex *vr, integer *ldvr, integer *mm, integer *m, complex *work, 
	real *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ctrexc)(char *compq, integer *n, complex *t, integer *
	ldt, complex *q, integer *ldq, integer *ifst, integer *ilst, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(ctrrfs)(char *uplo, char *trans, char *diag, integer *n, 
	integer *nrhs, complex *a, integer *lda, complex *b, integer *ldb, 
	complex *x, integer *ldx, real *ferr, real *berr, complex *work, real 
	*rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ctrsen)(char *job, char *compq, logical *select, integer 
	*n, complex *t, integer *ldt, complex *q, integer *ldq, complex *w, 
	integer *m, real *s, real *sep, complex *work, integer *lwork, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ctrsna)(char *job, char *howmny, logical *select, 
	integer *n, complex *t, integer *ldt, complex *vl, integer *ldvl, 
	complex *vr, integer *ldvr, real *s, real *sep, integer *mm, integer *
	m, complex *work, integer *ldwork, real *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ctrsyl)(char *trana, char *tranb, integer *isgn, integer 
	*m, integer *n, complex *a, integer *lda, complex *b, integer *ldb, 
	complex *c__, integer *ldc, real *scale, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ctrti2)(char *uplo, char *diag, integer *n, complex *a, 
	integer *lda, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ctrtri)(char *uplo, char *diag, integer *n, complex *a, 
	integer *lda, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ctrtrs)(char *uplo, char *trans, char *diag, integer *n, 
	integer *nrhs, complex *a, integer *lda, complex *b, integer *ldb, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ctzrqf)(integer *m, integer *n, complex *a, integer *lda,
	 complex *tau, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ctzrzf)(integer *m, integer *n, complex *a, integer *lda,
	 complex *tau, complex *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cung2l)(integer *m, integer *n, integer *k, complex *a, 
	integer *lda, complex *tau, complex *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cung2r)(integer *m, integer *n, integer *k, complex *a, 
	integer *lda, complex *tau, complex *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cungbr)(char *vect, integer *m, integer *n, integer *k, 
	complex *a, integer *lda, complex *tau, complex *work, integer *lwork,
	 integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cunghr)(integer *n, integer *ilo, integer *ihi, complex *
	a, integer *lda, complex *tau, complex *work, integer *lwork, integer 
	*info);
 
/* Subroutine */ int LAPACKMANGLE(cungl2)(integer *m, integer *n, integer *k, complex *a, 
	integer *lda, complex *tau, complex *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cunglq)(integer *m, integer *n, integer *k, complex *a, 
	integer *lda, complex *tau, complex *work, integer *lwork, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(cungql)(integer *m, integer *n, integer *k, complex *a, 
	integer *lda, complex *tau, complex *work, integer *lwork, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(cungqr)(integer *m, integer *n, integer *k, complex *a, 
	integer *lda, complex *tau, complex *work, integer *lwork, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(cungr2)(integer *m, integer *n, integer *k, complex *a, 
	integer *lda, complex *tau, complex *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cungrq)(integer *m, integer *n, integer *k, complex *a, 
	integer *lda, complex *tau, complex *work, integer *lwork, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(cungtr)(char *uplo, integer *n, complex *a, integer *lda,
	 complex *tau, complex *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cunm2l)(char *side, char *trans, integer *m, integer *n, 
	integer *k, complex *a, integer *lda, complex *tau, complex *c__, 
	integer *ldc, complex *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cunm2r)(char *side, char *trans, integer *m, integer *n, 
	integer *k, complex *a, integer *lda, complex *tau, complex *c__, 
	integer *ldc, complex *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cunmbr)(char *vect, char *side, char *trans, integer *m, 
	integer *n, integer *k, complex *a, integer *lda, complex *tau, 
	complex *c__, integer *ldc, complex *work, integer *lwork, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(cunmhr)(char *side, char *trans, integer *m, integer *n, 
	integer *ilo, integer *ihi, complex *a, integer *lda, complex *tau, 
	complex *c__, integer *ldc, complex *work, integer *lwork, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(cunml2)(char *side, char *trans, integer *m, integer *n, 
	integer *k, complex *a, integer *lda, complex *tau, complex *c__, 
	integer *ldc, complex *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cunmlq)(char *side, char *trans, integer *m, integer *n, 
	integer *k, complex *a, integer *lda, complex *tau, complex *c__, 
	integer *ldc, complex *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cunmql)(char *side, char *trans, integer *m, integer *n, 
	integer *k, complex *a, integer *lda, complex *tau, complex *c__, 
	integer *ldc, complex *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cunmqr)(char *side, char *trans, integer *m, integer *n, 
	integer *k, complex *a, integer *lda, complex *tau, complex *c__, 
	integer *ldc, complex *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cunmr2)(char *side, char *trans, integer *m, integer *n, 
	integer *k, complex *a, integer *lda, complex *tau, complex *c__, 
	integer *ldc, complex *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cunmr3)(char *side, char *trans, integer *m, integer *n, 
	integer *k, integer *l, complex *a, integer *lda, complex *tau, 
	complex *c__, integer *ldc, complex *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cunmrq)(char *side, char *trans, integer *m, integer *n, 
	integer *k, complex *a, integer *lda, complex *tau, complex *c__, 
	integer *ldc, complex *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cunmrz)(char *side, char *trans, integer *m, integer *n, 
	integer *k, integer *l, complex *a, integer *lda, complex *tau, 
	complex *c__, integer *ldc, complex *work, integer *lwork, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(cunmtr)(char *side, char *uplo, char *trans, integer *m, 
	integer *n, complex *a, integer *lda, complex *tau, complex *c__, 
	integer *ldc, complex *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cupgtr)(char *uplo, integer *n, complex *ap, complex *
	tau, complex *q, integer *ldq, complex *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(cupmtr)(char *side, char *uplo, char *trans, integer *m, 
	integer *n, complex *ap, complex *tau, complex *c__, integer *ldc, 
	complex *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dbdsdc)(char *uplo, char *compq, integer *n, doublereal *
	d__, doublereal *e, doublereal *u, integer *ldu, doublereal *vt, 
	integer *ldvt, doublereal *q, integer *iq, doublereal *work, integer *
	iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dbdsqr)(char *uplo, integer *n, integer *ncvt, integer *
	nru, integer *ncc, doublereal *d__, doublereal *e, doublereal *vt, 
	integer *ldvt, doublereal *u, integer *ldu, doublereal *c__, integer *
	ldc, doublereal *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ddisna)(char *job, integer *m, integer *n, doublereal *
	d__, doublereal *sep, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dgbbrd)(char *vect, integer *m, integer *n, integer *ncc,
	 integer *kl, integer *ku, doublereal *ab, integer *ldab, doublereal *
	d__, doublereal *e, doublereal *q, integer *ldq, doublereal *pt, 
	integer *ldpt, doublereal *c__, integer *ldc, doublereal *work, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dgbcon)(char *norm, integer *n, integer *kl, integer *ku,
	 doublereal *ab, integer *ldab, integer *ipiv, doublereal *anorm, 
	doublereal *rcond, doublereal *work, integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dgbequ)(integer *m, integer *n, integer *kl, integer *ku,
	 doublereal *ab, integer *ldab, doublereal *r__, doublereal *c__, 
	doublereal *rowcnd, doublereal *colcnd, doublereal *amax, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(dgbrfs)(char *trans, integer *n, integer *kl, integer *
	ku, integer *nrhs, doublereal *ab, integer *ldab, doublereal *afb, 
	integer *ldafb, integer *ipiv, doublereal *b, integer *ldb, 
	doublereal *x, integer *ldx, doublereal *ferr, doublereal *berr, 
	doublereal *work, integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dgbsv)(integer *n, integer *kl, integer *ku, integer *
	nrhs, doublereal *ab, integer *ldab, integer *ipiv, doublereal *b, 
	integer *ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dgbsvx)(char *fact, char *trans, integer *n, integer *kl,
	 integer *ku, integer *nrhs, doublereal *ab, integer *ldab, 
	doublereal *afb, integer *ldafb, integer *ipiv, char *equed, 
	doublereal *r__, doublereal *c__, doublereal *b, integer *ldb, 
	doublereal *x, integer *ldx, doublereal *rcond, doublereal *ferr, 
	doublereal *berr, doublereal *work, integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dgbtf2)(integer *m, integer *n, integer *kl, integer *ku,
	 doublereal *ab, integer *ldab, integer *ipiv, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dgbtrf)(integer *m, integer *n, integer *kl, integer *ku,
	 doublereal *ab, integer *ldab, integer *ipiv, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dgbtrs)(char *trans, integer *n, integer *kl, integer *
	ku, integer *nrhs, doublereal *ab, integer *ldab, integer *ipiv, 
	doublereal *b, integer *ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dgebak)(char *job, char *side, integer *n, integer *ilo, 
	integer *ihi, doublereal *scale, integer *m, doublereal *v, integer *
	ldv, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dgebal)(char *job, integer *n, doublereal *a, integer *
	lda, integer *ilo, integer *ihi, doublereal *scale, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dgebd2)(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *d__, doublereal *e, doublereal *tauq, doublereal *
	taup, doublereal *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dgebrd)(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *d__, doublereal *e, doublereal *tauq, doublereal *
	taup, doublereal *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dgecon)(char *norm, integer *n, doublereal *a, integer *
	lda, doublereal *anorm, doublereal *rcond, doublereal *work, integer *
	iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dgeequ)(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *r__, doublereal *c__, doublereal *rowcnd, doublereal 
	*colcnd, doublereal *amax, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dgees)(char *jobvs, char *sort, L_fp select, integer *n, 
	doublereal *a, integer *lda, integer *sdim, doublereal *wr, 
	doublereal *wi, doublereal *vs, integer *ldvs, doublereal *work, 
	integer *lwork, logical *bwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dgeesx)(char *jobvs, char *sort, L_fp select, char *
	sense, integer *n, doublereal *a, integer *lda, integer *sdim, 
	doublereal *wr, doublereal *wi, doublereal *vs, integer *ldvs, 
	doublereal *rconde, doublereal *rcondv, doublereal *work, integer *
	lwork, integer *iwork, integer *liwork, logical *bwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dgeev)(char *jobvl, char *jobvr, integer *n, doublereal *
	a, integer *lda, doublereal *wr, doublereal *wi, doublereal *vl, 
	integer *ldvl, doublereal *vr, integer *ldvr, doublereal *work, 
	integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dgeevx)(char *balanc, char *jobvl, char *jobvr, char *
	sense, integer *n, doublereal *a, integer *lda, doublereal *wr, 
	doublereal *wi, doublereal *vl, integer *ldvl, doublereal *vr, 
	integer *ldvr, integer *ilo, integer *ihi, doublereal *scale, 
	doublereal *abnrm, doublereal *rconde, doublereal *rcondv, doublereal 
	*work, integer *lwork, integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dgegs)(char *jobvsl, char *jobvsr, integer *n, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	alphar, doublereal *alphai, doublereal *beta, doublereal *vsl, 
	integer *ldvsl, doublereal *vsr, integer *ldvsr, doublereal *work, 
	integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dgegv)(char *jobvl, char *jobvr, integer *n, doublereal *
	a, integer *lda, doublereal *b, integer *ldb, doublereal *alphar, 
	doublereal *alphai, doublereal *beta, doublereal *vl, integer *ldvl, 
	doublereal *vr, integer *ldvr, doublereal *work, integer *lwork, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dgehd2)(integer *n, integer *ilo, integer *ihi, 
	doublereal *a, integer *lda, doublereal *tau, doublereal *work, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dgehrd)(integer *n, integer *ilo, integer *ihi, 
	doublereal *a, integer *lda, doublereal *tau, doublereal *work, 
	integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dgelq2)(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *tau, doublereal *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dgelqf)(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *tau, doublereal *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dgels)(char *trans, integer *m, integer *n, integer *
	nrhs, doublereal *a, integer *lda, doublereal *b, integer *ldb, 
	doublereal *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dgelsd)(integer *m, integer *n, integer *nrhs, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	s, doublereal *rcond, integer *rank, doublereal *work, integer *lwork,
	 integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dgelss)(integer *m, integer *n, integer *nrhs, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	s, doublereal *rcond, integer *rank, doublereal *work, integer *lwork,
	 integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dgelsx)(integer *m, integer *n, integer *nrhs, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, integer *
	jpvt, doublereal *rcond, integer *rank, doublereal *work, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(dgelsy)(integer *m, integer *n, integer *nrhs, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, integer *
	jpvt, doublereal *rcond, integer *rank, doublereal *work, integer *
	lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dgeql2)(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *tau, doublereal *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dgeqlf)(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *tau, doublereal *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dgeqp3)(integer *m, integer *n, doublereal *a, integer *
	lda, integer *jpvt, doublereal *tau, doublereal *work, integer *lwork,
	 integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dgeqpf)(integer *m, integer *n, doublereal *a, integer *
	lda, integer *jpvt, doublereal *tau, doublereal *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dgeqr2)(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *tau, doublereal *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dgeqrf)(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *tau, doublereal *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dgerfs)(char *trans, integer *n, integer *nrhs, 
	doublereal *a, integer *lda, doublereal *af, integer *ldaf, integer *
	ipiv, doublereal *b, integer *ldb, doublereal *x, integer *ldx, 
	doublereal *ferr, doublereal *berr, doublereal *work, integer *iwork, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dgerq2)(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *tau, doublereal *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dgerqf)(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *tau, doublereal *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dgesc2)(integer *n, doublereal *a, integer *lda, 
	doublereal *rhs, integer *ipiv, integer *jpiv, doublereal *scale);
 
/* Subroutine */ int LAPACKMANGLE(dgesdd)(char *jobz, integer *m, integer *n, doublereal *
	a, integer *lda, doublereal *s, doublereal *u, integer *ldu, 
	doublereal *vt, integer *ldvt, doublereal *work, integer *lwork, 
	integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dgesv)(integer *n, integer *nrhs, doublereal *a, integer 
	*lda, integer *ipiv, doublereal *b, integer *ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dgesvd)(char *jobu, char *jobvt, integer *m, integer *n, 
	doublereal *a, integer *lda, doublereal *s, doublereal *u, integer *
	ldu, doublereal *vt, integer *ldvt, doublereal *work, integer *lwork, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dgesvx)(char *fact, char *trans, integer *n, integer *
	nrhs, doublereal *a, integer *lda, doublereal *af, integer *ldaf, 
	integer *ipiv, char *equed, doublereal *r__, doublereal *c__, 
	doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *
	rcond, doublereal *ferr, doublereal *berr, doublereal *work, integer *
	iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dgetc2)(integer *n, doublereal *a, integer *lda, integer 
	*ipiv, integer *jpiv, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dgetf2)(integer *m, integer *n, doublereal *a, integer *
	lda, integer *ipiv, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dgetrf)(integer *m, integer *n, doublereal *a, integer *
	lda, integer *ipiv, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dgetri)(integer *n, doublereal *a, integer *lda, integer 
	*ipiv, doublereal *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dgetrs)(char *trans, integer *n, integer *nrhs, 
	doublereal *a, integer *lda, integer *ipiv, doublereal *b, integer *
	ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dggbak)(char *job, char *side, integer *n, integer *ilo, 
	integer *ihi, doublereal *lscale, doublereal *rscale, integer *m, 
	doublereal *v, integer *ldv, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dggbal)(char *job, integer *n, doublereal *a, integer *
	lda, doublereal *b, integer *ldb, integer *ilo, integer *ihi, 
	doublereal *lscale, doublereal *rscale, doublereal *work, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(dgges)(char *jobvsl, char *jobvsr, char *sort, L_fp 
	delctg, integer *n, doublereal *a, integer *lda, doublereal *b, 
	integer *ldb, integer *sdim, doublereal *alphar, doublereal *alphai, 
	doublereal *beta, doublereal *vsl, integer *ldvsl, doublereal *vsr, 
	integer *ldvsr, doublereal *work, integer *lwork, logical *bwork, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dggesx)(char *jobvsl, char *jobvsr, char *sort, L_fp 
	delctg, char *sense, integer *n, doublereal *a, integer *lda, 
	doublereal *b, integer *ldb, integer *sdim, doublereal *alphar, 
	doublereal *alphai, doublereal *beta, doublereal *vsl, integer *ldvsl,
	 doublereal *vsr, integer *ldvsr, doublereal *rconde, doublereal *
	rcondv, doublereal *work, integer *lwork, integer *iwork, integer *
	liwork, logical *bwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dggev)(char *jobvl, char *jobvr, integer *n, doublereal *
	a, integer *lda, doublereal *b, integer *ldb, doublereal *alphar, 
	doublereal *alphai, doublereal *beta, doublereal *vl, integer *ldvl, 
	doublereal *vr, integer *ldvr, doublereal *work, integer *lwork, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dggevx)(char *balanc, char *jobvl, char *jobvr, char *
	sense, integer *n, doublereal *a, integer *lda, doublereal *b, 
	integer *ldb, doublereal *alphar, doublereal *alphai, doublereal *
	beta, doublereal *vl, integer *ldvl, doublereal *vr, integer *ldvr, 
	integer *ilo, integer *ihi, doublereal *lscale, doublereal *rscale, 
	doublereal *abnrm, doublereal *bbnrm, doublereal *rconde, doublereal *
	rcondv, doublereal *work, integer *lwork, integer *iwork, logical *
	bwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dggglm)(integer *n, integer *m, integer *p, doublereal *
	a, integer *lda, doublereal *b, integer *ldb, doublereal *d__, 
	doublereal *x, doublereal *y, doublereal *work, integer *lwork, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dgghrd)(char *compq, char *compz, integer *n, integer *
	ilo, integer *ihi, doublereal *a, integer *lda, doublereal *b, 
	integer *ldb, doublereal *q, integer *ldq, doublereal *z__, integer *
	ldz, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dgglse)(integer *m, integer *n, integer *p, doublereal *
	a, integer *lda, doublereal *b, integer *ldb, doublereal *c__, 
	doublereal *d__, doublereal *x, doublereal *work, integer *lwork, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dggqrf)(integer *n, integer *m, integer *p, doublereal *
	a, integer *lda, doublereal *taua, doublereal *b, integer *ldb, 
	doublereal *taub, doublereal *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dggrqf)(integer *m, integer *p, integer *n, doublereal *
	a, integer *lda, doublereal *taua, doublereal *b, integer *ldb, 
	doublereal *taub, doublereal *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dggsvd)(char *jobu, char *jobv, char *jobq, integer *m, 
	integer *n, integer *p, integer *k, integer *l, doublereal *a, 
	integer *lda, doublereal *b, integer *ldb, doublereal *alpha, 
	doublereal *beta, doublereal *u, integer *ldu, doublereal *v, integer 
	*ldv, doublereal *q, integer *ldq, doublereal *work, integer *iwork, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dggsvp)(char *jobu, char *jobv, char *jobq, integer *m, 
	integer *p, integer *n, doublereal *a, integer *lda, doublereal *b, 
	integer *ldb, doublereal *tola, doublereal *tolb, integer *k, integer 
	*l, doublereal *u, integer *ldu, doublereal *v, integer *ldv, 
	doublereal *q, integer *ldq, integer *iwork, doublereal *tau, 
	doublereal *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dgtcon)(char *norm, integer *n, doublereal *dl, 
	doublereal *d__, doublereal *du, doublereal *du2, integer *ipiv, 
	doublereal *anorm, doublereal *rcond, doublereal *work, integer *
	iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dgtrfs)(char *trans, integer *n, integer *nrhs, 
	doublereal *dl, doublereal *d__, doublereal *du, doublereal *dlf, 
	doublereal *df, doublereal *duf, doublereal *du2, integer *ipiv, 
	doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *
	ferr, doublereal *berr, doublereal *work, integer *iwork, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(dgtsv)(integer *n, integer *nrhs, doublereal *dl, 
	doublereal *d__, doublereal *du, doublereal *b, integer *ldb, integer 
	*info);
 
/* Subroutine */ int LAPACKMANGLE(dgtsvx)(char *fact, char *trans, integer *n, integer *
	nrhs, doublereal *dl, doublereal *d__, doublereal *du, doublereal *
	dlf, doublereal *df, doublereal *duf, doublereal *du2, integer *ipiv, 
	doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *
	rcond, doublereal *ferr, doublereal *berr, doublereal *work, integer *
	iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dgttrf)(integer *n, doublereal *dl, doublereal *d__, 
	doublereal *du, doublereal *du2, integer *ipiv, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dgttrs)(char *trans, integer *n, integer *nrhs, 
	doublereal *dl, doublereal *d__, doublereal *du, doublereal *du2, 
	integer *ipiv, doublereal *b, integer *ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dgtts2)(integer *itrans, integer *n, integer *nrhs, 
	doublereal *dl, doublereal *d__, doublereal *du, doublereal *du2, 
	integer *ipiv, doublereal *b, integer *ldb);
 
/* Subroutine */ int LAPACKMANGLE(dhgeqz)(char *job, char *compq, char *compz, integer *n, 
	integer *ilo, integer *ihi, doublereal *a, integer *lda, doublereal *
	b, integer *ldb, doublereal *alphar, doublereal *alphai, doublereal *
	beta, doublereal *q, integer *ldq, doublereal *z__, integer *ldz, 
	doublereal *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dhsein)(char *side, char *eigsrc, char *initv, logical *
	select, integer *n, doublereal *h__, integer *ldh, doublereal *wr, 
	doublereal *wi, doublereal *vl, integer *ldvl, doublereal *vr, 
	integer *ldvr, integer *mm, integer *m, doublereal *work, integer *
	ifaill, integer *ifailr, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dhseqr)(char *job, char *compz, integer *n, integer *ilo,
	 integer *ihi, doublereal *h__, integer *ldh, doublereal *wr, 
	doublereal *wi, doublereal *z__, integer *ldz, doublereal *work, 
	integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dlabad)(doublereal *small, doublereal *large);
 
/* Subroutine */ int LAPACKMANGLE(dlabrd)(integer *m, integer *n, integer *nb, doublereal *
	a, integer *lda, doublereal *d__, doublereal *e, doublereal *tauq, 
	doublereal *taup, doublereal *x, integer *ldx, doublereal *y, integer 
	*ldy);
 
/* Subroutine */ int LAPACKMANGLE(dlacon)(integer *n, doublereal *v, doublereal *x, 
	integer *isgn, doublereal *est, integer *kase);
 
/* Subroutine */ int LAPACKMANGLE(dlacpy)(char *uplo, integer *m, integer *n, doublereal *
	a, integer *lda, doublereal *b, integer *ldb);
 
/* Subroutine */ int LAPACKMANGLE(dladiv)(doublereal *a, doublereal *b, doublereal *c__, 
	doublereal *d__, doublereal *p, doublereal *q);
 
/* Subroutine */ int LAPACKMANGLE(dlae2)(doublereal *a, doublereal *b, doublereal *c__, 
	doublereal *rt1, doublereal *rt2);
 
/* Subroutine */ int LAPACKMANGLE(dlaebz)(integer *ijob, integer *nitmax, integer *n, 
	integer *mmax, integer *minp, integer *nbmin, doublereal *abstol, 
	doublereal *reltol, doublereal *pivmin, doublereal *d__, doublereal *
	e, doublereal *e2, integer *nval, doublereal *ab, doublereal *c__, 
	integer *mout, integer *nab, doublereal *work, integer *iwork, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dlaed0)(integer *icompq, integer *qsiz, integer *n, 
	doublereal *d__, doublereal *e, doublereal *q, integer *ldq, 
	doublereal *qstore, integer *ldqs, doublereal *work, integer *iwork, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dlaed1)(integer *n, doublereal *d__, doublereal *q, 
	integer *ldq, integer *indxq, doublereal *rho, integer *cutpnt, 
	doublereal *work, integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dlaed2)(integer *k, integer *n, integer *n1, doublereal *
	d__, doublereal *q, integer *ldq, integer *indxq, doublereal *rho, 
	doublereal *z__, doublereal *dlamda, doublereal *w, doublereal *q2, 
	integer *indx, integer *indxc, integer *indxp, integer *coltyp, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dlaed3)(integer *k, integer *n, integer *n1, doublereal *
	d__, doublereal *q, integer *ldq, doublereal *rho, doublereal *dlamda,
	 doublereal *q2, integer *indx, integer *ctot, doublereal *w, 
	doublereal *s, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dlaed4)(integer *n, integer *i__, doublereal *d__, 
	doublereal *z__, doublereal *delta, doublereal *rho, doublereal *dlam,
	 integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dlaed5)(integer *i__, doublereal *d__, doublereal *z__, 
	doublereal *delta, doublereal *rho, doublereal *dlam);
 
/* Subroutine */ int LAPACKMANGLE(dlaed6)(integer *kniter, logical *orgati, doublereal *
	rho, doublereal *d__, doublereal *z__, doublereal *finit, doublereal *
	tau, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dlaed7)(integer *icompq, integer *n, integer *qsiz, 
	integer *tlvls, integer *curlvl, integer *curpbm, doublereal *d__, 
	doublereal *q, integer *ldq, integer *indxq, doublereal *rho, integer 
	*cutpnt, doublereal *qstore, integer *qptr, integer *prmptr, integer *
	perm, integer *givptr, integer *givcol, doublereal *givnum, 
	doublereal *work, integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dlaed8)(integer *icompq, integer *k, integer *n, integer 
	*qsiz, doublereal *d__, doublereal *q, integer *ldq, integer *indxq, 
	doublereal *rho, integer *cutpnt, doublereal *z__, doublereal *dlamda,
	 doublereal *q2, integer *ldq2, doublereal *w, integer *perm, integer 
	*givptr, integer *givcol, doublereal *givnum, integer *indxp, integer 
	*indx, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dlaed9)(integer *k, integer *kstart, integer *kstop, 
	integer *n, doublereal *d__, doublereal *q, integer *ldq, doublereal *
	rho, doublereal *dlamda, doublereal *w, doublereal *s, integer *lds, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dlaeda)(integer *n, integer *tlvls, integer *curlvl, 
	integer *curpbm, integer *prmptr, integer *perm, integer *givptr, 
	integer *givcol, doublereal *givnum, doublereal *q, integer *qptr, 
	doublereal *z__, doublereal *ztemp, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dlaein)(logical *rightv, logical *noinit, integer *n, 
	doublereal *h__, integer *ldh, doublereal *wr, doublereal *wi, 
	doublereal *vr, doublereal *vi, doublereal *b, integer *ldb, 
	doublereal *work, doublereal *eps3, doublereal *smlnum, doublereal *
	bignum, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dlaev2)(doublereal *a, doublereal *b, doublereal *c__, 
	doublereal *rt1, doublereal *rt2, doublereal *cs1, doublereal *sn1);
 
/* Subroutine */ int LAPACKMANGLE(dlaexc)(logical *wantq, integer *n, doublereal *t, 
	integer *ldt, doublereal *q, integer *ldq, integer *j1, integer *n1, 
	integer *n2, doublereal *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dlag2)(doublereal *a, integer *lda, doublereal *b, 
	integer *ldb, doublereal *safmin, doublereal *scale1, doublereal *
	scale2, doublereal *wr1, doublereal *wr2, doublereal *wi);
 
/* Subroutine */ int LAPACKMANGLE(dlags2)(logical *upper, doublereal *a1, doublereal *a2, 
	doublereal *a3, doublereal *b1, doublereal *b2, doublereal *b3, 
	doublereal *csu, doublereal *snu, doublereal *csv, doublereal *snv, 
	doublereal *csq, doublereal *snq);
 
/* Subroutine */ int LAPACKMANGLE(dlagtf)(integer *n, doublereal *a, doublereal *lambda, 
	doublereal *b, doublereal *c__, doublereal *tol, doublereal *d__, 
	integer *in, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dlagtm)(char *trans, integer *n, integer *nrhs, 
	doublereal *alpha, doublereal *dl, doublereal *d__, doublereal *du, 
	doublereal *x, integer *ldx, doublereal *beta, doublereal *b, integer 
	*ldb);
 
/* Subroutine */ int LAPACKMANGLE(dlagts)(integer *job, integer *n, doublereal *a, 
	doublereal *b, doublereal *c__, doublereal *d__, integer *in, 
	doublereal *y, doublereal *tol, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dlagv2)(doublereal *a, integer *lda, doublereal *b, 
	integer *ldb, doublereal *alphar, doublereal *alphai, doublereal *
	beta, doublereal *csl, doublereal *snl, doublereal *csr, doublereal *
	snr);
 
/* Subroutine */ int LAPACKMANGLE(dlahqr)(logical *wantt, logical *wantz, integer *n, 
	integer *ilo, integer *ihi, doublereal *h__, integer *ldh, doublereal 
	*wr, doublereal *wi, integer *iloz, integer *ihiz, doublereal *z__, 
	integer *ldz, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dlahrd)(integer *n, integer *k, integer *nb, doublereal *
	a, integer *lda, doublereal *tau, doublereal *t, integer *ldt, 
	doublereal *y, integer *ldy);
 
/* Subroutine */ int LAPACKMANGLE(dlaic1)(integer *job, integer *j, doublereal *x, 
	doublereal *sest, doublereal *w, doublereal *gamma, doublereal *
	sestpr, doublereal *s, doublereal *c__);
 
/* Subroutine */ int LAPACKMANGLE(dlaln2)(logical *ltrans, integer *na, integer *nw, 
	doublereal *smin, doublereal *ca, doublereal *a, integer *lda, 
	doublereal *d1, doublereal *d2, doublereal *b, integer *ldb, 
	doublereal *wr, doublereal *wi, doublereal *x, integer *ldx, 
	doublereal *scale, doublereal *xnorm, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dlals0)(integer *icompq, integer *nl, integer *nr, 
	integer *sqre, integer *nrhs, doublereal *b, integer *ldb, doublereal 
	*bx, integer *ldbx, integer *perm, integer *givptr, integer *givcol, 
	integer *ldgcol, doublereal *givnum, integer *ldgnum, doublereal *
	poles, doublereal *difl, doublereal *difr, doublereal *z__, integer *
	k, doublereal *c__, doublereal *s, doublereal *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dlalsa)(integer *icompq, integer *smlsiz, integer *n, 
	integer *nrhs, doublereal *b, integer *ldb, doublereal *bx, integer *
	ldbx, doublereal *u, integer *ldu, doublereal *vt, integer *k, 
	doublereal *difl, doublereal *difr, doublereal *z__, doublereal *
	poles, integer *givptr, integer *givcol, integer *ldgcol, integer *
	perm, doublereal *givnum, doublereal *c__, doublereal *s, doublereal *
	work, integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dlalsd)(char *uplo, integer *smlsiz, integer *n, integer 
	*nrhs, doublereal *d__, doublereal *e, doublereal *b, integer *ldb, 
	doublereal *rcond, integer *rank, doublereal *work, integer *iwork, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dlamc1)(integer *beta, integer *t, logical *rnd, logical 
	*ieee1);
 
/* Subroutine */ int LAPACKMANGLE(dlamc2)(integer *beta, integer *t, logical *rnd, 
	doublereal *eps, integer *emin, doublereal *rmin, integer *emax, 
	doublereal *rmax);
 
/* Subroutine */ int LAPACKMANGLE(dlamc4)(integer *emin, doublereal *start, integer *base);
 
/* Subroutine */ int LAPACKMANGLE(dlamc5)(integer *beta, integer *p, integer *emin, 
	logical *ieee, integer *emax, doublereal *rmax);
 
/* Subroutine */ int LAPACKMANGLE(dlamrg)(integer *n1, integer *n2, doublereal *a, integer 
	*dtrd1, integer *dtrd2, integer *index);
 
/* Subroutine */ int LAPACKMANGLE(dlanv2)(doublereal *a, doublereal *b, doublereal *c__, 
	doublereal *d__, doublereal *rt1r, doublereal *rt1i, doublereal *rt2r,
	 doublereal *rt2i, doublereal *cs, doublereal *sn);
 
/* Subroutine */ int LAPACKMANGLE(dlapll)(integer *n, doublereal *x, integer *incx, 
	doublereal *y, integer *incy, doublereal *ssmin);
 
/* Subroutine */ int LAPACKMANGLE(dlapmt)(logical *forwrd, integer *m, integer *n, 
	doublereal *x, integer *ldx, integer *k);
 
/* Subroutine */ int LAPACKMANGLE(dlaqgb)(integer *m, integer *n, integer *kl, integer *ku,
	 doublereal *ab, integer *ldab, doublereal *r__, doublereal *c__, 
	doublereal *rowcnd, doublereal *colcnd, doublereal *amax, char *equed);
 
/* Subroutine */ int LAPACKMANGLE(dlaqge)(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *r__, doublereal *c__, doublereal *rowcnd, doublereal 
	*colcnd, doublereal *amax, char *equed);
 
/* Subroutine */ int LAPACKMANGLE(dlaqp2)(integer *m, integer *n, integer *offset, 
	doublereal *a, integer *lda, integer *jpvt, doublereal *tau, 
	doublereal *vn1, doublereal *vn2, doublereal *work);
 
/* Subroutine */ int LAPACKMANGLE(dlaqps)(integer *m, integer *n, integer *offset, integer 
	*nb, integer *kb, doublereal *a, integer *lda, integer *jpvt, 
	doublereal *tau, doublereal *vn1, doublereal *vn2, doublereal *auxv, 
	doublereal *f, integer *ldf);
 
/* Subroutine */ int LAPACKMANGLE(dlaqsb)(char *uplo, integer *n, integer *kd, doublereal *
	ab, integer *ldab, doublereal *s, doublereal *scond, doublereal *amax,
	 char *equed);
 
/* Subroutine */ int LAPACKMANGLE(dlaqsp)(char *uplo, integer *n, doublereal *ap, 
	doublereal *s, doublereal *scond, doublereal *amax, char *equed);
 
/* Subroutine */ int LAPACKMANGLE(dlaqsy)(char *uplo, integer *n, doublereal *a, integer *
	lda, doublereal *s, doublereal *scond, doublereal *amax, char *equed);
 
/* Subroutine */ int LAPACKMANGLE(dlaqtr)(logical *ltran, logical *lreal, integer *n, 
	doublereal *t, integer *ldt, doublereal *b, doublereal *w, doublereal 
	*scale, doublereal *x, doublereal *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dlar1v)(integer *n, integer *b1, integer *bn, doublereal 
	*sigma, doublereal *d__, doublereal *l, doublereal *ld, doublereal *
	lld, doublereal *gersch, doublereal *z__, doublereal *ztz, doublereal 
	*mingma, integer *r__, integer *isuppz, doublereal *work);
 
/* Subroutine */ int LAPACKMANGLE(dlar2v)(integer *n, doublereal *x, doublereal *y, 
	doublereal *z__, integer *incx, doublereal *c__, doublereal *s, 
	integer *incc);
 
/* Subroutine */ int LAPACKMANGLE(dlarf)(char *side, integer *m, integer *n, doublereal *v,
	 integer *incv, doublereal *tau, doublereal *c__, integer *ldc, 
	doublereal *work);
 
/* Subroutine */ int LAPACKMANGLE(dlarfb)(char *side, char *trans, char *direct, char *
	storev, integer *m, integer *n, integer *k, doublereal *v, integer *
	ldv, doublereal *t, integer *ldt, doublereal *c__, integer *ldc, 
	doublereal *work, integer *ldwork);
 
/* Subroutine */ int LAPACKMANGLE(dlarfg)(integer *n, doublereal *alpha, doublereal *x, 
	integer *incx, doublereal *tau);
 
/* Subroutine */ int LAPACKMANGLE(dlarft)(char *direct, char *storev, integer *n, integer *
	k, doublereal *v, integer *ldv, doublereal *tau, doublereal *t, 
	integer *ldt);
 
/* Subroutine */ int LAPACKMANGLE(dlarfx)(char *side, integer *m, integer *n, doublereal *
	v, doublereal *tau, doublereal *c__, integer *ldc, doublereal *work);
 
/* Subroutine */ int LAPACKMANGLE(dlargv)(integer *n, doublereal *x, integer *incx, 
	doublereal *y, integer *incy, doublereal *c__, integer *incc);
 
/* Subroutine */ int LAPACKMANGLE(dlarnv)(integer *idist, integer *iseed, integer *n, 
	doublereal *x);
 
/* Subroutine */ int LAPACKMANGLE(dlarrb)(integer *n, doublereal *d__, doublereal *l, 
	doublereal *ld, doublereal *lld, integer *ifirst, integer *ilast, 
	doublereal *sigma, doublereal *reltol, doublereal *w, doublereal *
	wgap, doublereal *werr, doublereal *work, integer *iwork, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(dlarre)(integer *n, doublereal *d__, doublereal *e, 
	doublereal *tol, integer *nsplit, integer *isplit, integer *m, 
	doublereal *w, doublereal *woff, doublereal *gersch, doublereal *work,
	 integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dlarrf)(integer *n, doublereal *d__, doublereal *l, 
	doublereal *ld, doublereal *lld, integer *ifirst, integer *ilast, 
	doublereal *w, doublereal *dplus, doublereal *lplus, doublereal *work,
	 integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dlarrv)(integer *n, doublereal *d__, doublereal *l, 
	integer *isplit, integer *m, doublereal *w, integer *iblock, 
	doublereal *gersch, doublereal *tol, doublereal *z__, integer *ldz, 
	integer *isuppz, doublereal *work, integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dlartg)(doublereal *f, doublereal *g, doublereal *cs, 
	doublereal *sn, doublereal *r__);
 
/* Subroutine */ int LAPACKMANGLE(dlartv)(integer *n, doublereal *x, integer *incx, 
	doublereal *y, integer *incy, doublereal *c__, doublereal *s, integer 
	*incc);
 
/* Subroutine */ int LAPACKMANGLE(dlaruv)(integer *iseed, integer *n, doublereal *x);
 
/* Subroutine */ int LAPACKMANGLE(dlarz)(char *side, integer *m, integer *n, integer *l, 
	doublereal *v, integer *incv, doublereal *tau, doublereal *c__, 
	integer *ldc, doublereal *work);
 
/* Subroutine */ int LAPACKMANGLE(dlarzb)(char *side, char *trans, char *direct, char *
	storev, integer *m, integer *n, integer *k, integer *l, doublereal *v,
	 integer *ldv, doublereal *t, integer *ldt, doublereal *c__, integer *
	ldc, doublereal *work, integer *ldwork);
 
/* Subroutine */ int LAPACKMANGLE(dlarzt)(char *direct, char *storev, integer *n, integer *
	k, doublereal *v, integer *ldv, doublereal *tau, doublereal *t, 
	integer *ldt);
 
/* Subroutine */ int LAPACKMANGLE(dlas2)(doublereal *f, doublereal *g, doublereal *h__, 
	doublereal *ssmin, doublereal *ssmax);
 
/* Subroutine */ int LAPACKMANGLE(dlascl)(char *type__, integer *kl, integer *ku, 
	doublereal *cfrom, doublereal *cto, integer *m, integer *n, 
	doublereal *a, integer *lda, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dlasd0)(integer *n, integer *sqre, doublereal *d__, 
	doublereal *e, doublereal *u, integer *ldu, doublereal *vt, integer *
	ldvt, integer *smlsiz, integer *iwork, doublereal *work, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(dlasd1)(integer *nl, integer *nr, integer *sqre, 
	doublereal *d__, doublereal *alpha, doublereal *beta, doublereal *u, 
	integer *ldu, doublereal *vt, integer *ldvt, integer *idxq, integer *
	iwork, doublereal *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dlasd2)(integer *nl, integer *nr, integer *sqre, integer 
	*k, doublereal *d__, doublereal *z__, doublereal *alpha, doublereal *
	beta, doublereal *u, integer *ldu, doublereal *vt, integer *ldvt, 
	doublereal *dsigma, doublereal *u2, integer *ldu2, doublereal *vt2, 
	integer *ldvt2, integer *idxp, integer *idx, integer *idxc, integer *
	idxq, integer *coltyp, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dlasd3)(integer *nl, integer *nr, integer *sqre, integer 
	*k, doublereal *d__, doublereal *q, integer *ldq, doublereal *dsigma, 
	doublereal *u, integer *ldu, doublereal *u2, integer *ldu2, 
	doublereal *vt, integer *ldvt, doublereal *vt2, integer *ldvt2, 
	integer *idxc, integer *ctot, doublereal *z__, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dlasd4)(integer *n, integer *i__, doublereal *d__, 
	doublereal *z__, doublereal *delta, doublereal *rho, doublereal *
	sigma, doublereal *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dlasd5)(integer *i__, doublereal *d__, doublereal *z__, 
	doublereal *delta, doublereal *rho, doublereal *dsigma, doublereal *
	work);
 
/* Subroutine */ int LAPACKMANGLE(dlasd6)(integer *icompq, integer *nl, integer *nr, 
	integer *sqre, doublereal *d__, doublereal *vf, doublereal *vl, 
	doublereal *alpha, doublereal *beta, integer *idxq, integer *perm, 
	integer *givptr, integer *givcol, integer *ldgcol, doublereal *givnum,
	 integer *ldgnum, doublereal *poles, doublereal *difl, doublereal *
	difr, doublereal *z__, integer *k, doublereal *c__, doublereal *s, 
	doublereal *work, integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dlasd7)(integer *icompq, integer *nl, integer *nr, 
	integer *sqre, integer *k, doublereal *d__, doublereal *z__, 
	doublereal *zw, doublereal *vf, doublereal *vfw, doublereal *vl, 
	doublereal *vlw, doublereal *alpha, doublereal *beta, doublereal *
	dsigma, integer *idx, integer *idxp, integer *idxq, integer *perm, 
	integer *givptr, integer *givcol, integer *ldgcol, doublereal *givnum,
	 integer *ldgnum, doublereal *c__, doublereal *s, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dlasd8)(integer *icompq, integer *k, doublereal *d__, 
	doublereal *z__, doublereal *vf, doublereal *vl, doublereal *difl, 
	doublereal *difr, integer *lddifr, doublereal *dsigma, doublereal *
	work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dlasd9)(integer *icompq, integer *ldu, integer *k, 
	doublereal *d__, doublereal *z__, doublereal *vf, doublereal *vl, 
	doublereal *difl, doublereal *difr, doublereal *dsigma, doublereal *
	work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dlasda)(integer *icompq, integer *smlsiz, integer *n, 
	integer *sqre, doublereal *d__, doublereal *e, doublereal *u, integer 
	*ldu, doublereal *vt, integer *k, doublereal *difl, doublereal *difr, 
	doublereal *z__, doublereal *poles, integer *givptr, integer *givcol, 
	integer *ldgcol, integer *perm, doublereal *givnum, doublereal *c__, 
	doublereal *s, doublereal *work, integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dlasdq)(char *uplo, integer *sqre, integer *n, integer *
	ncvt, integer *nru, integer *ncc, doublereal *d__, doublereal *e, 
	doublereal *vt, integer *ldvt, doublereal *u, integer *ldu, 
	doublereal *c__, integer *ldc, doublereal *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dlasdt)(integer *n, integer *lvl, integer *nd, integer *
	inode, integer *ndiml, integer *ndimr, integer *msub);
 
/* Subroutine */ int LAPACKMANGLE(dlaset)(char *uplo, integer *m, integer *n, doublereal *
	alpha, doublereal *beta, doublereal *a, integer *lda);
 
/* Subroutine */ int LAPACKMANGLE(dlasq1)(integer *n, doublereal *d__, doublereal *e, 
	doublereal *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dlasq2)(integer *n, doublereal *z__, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dlasq3)(integer *i0, integer *n0, doublereal *z__, 
	integer *pp, doublereal *dmin__, doublereal *sigma, doublereal *desig,
	 doublereal *qmax, integer *nfail, integer *iter, integer *ndiv, 
	logical *ieee);
 
/* Subroutine */ int LAPACKMANGLE(dlasq4)(integer *i0, integer *n0, doublereal *z__, 
	integer *pp, integer *n0in, doublereal *dmin__, doublereal *dmin1, 
	doublereal *dmin2, doublereal *dn, doublereal *dn1, doublereal *dn2, 
	doublereal *tau, integer *ttype);
 
/* Subroutine */ int LAPACKMANGLE(dlasq5)(integer *i0, integer *n0, doublereal *z__, 
	integer *pp, doublereal *tau, doublereal *dmin__, doublereal *dmin1, 
	doublereal *dmin2, doublereal *dn, doublereal *dnm1, doublereal *dnm2,
	 logical *ieee);
 
/* Subroutine */ int LAPACKMANGLE(dlasq6)(integer *i0, integer *n0, doublereal *z__, 
	integer *pp, doublereal *dmin__, doublereal *dmin1, doublereal *dmin2,
	 doublereal *dn, doublereal *dnm1, doublereal *dnm2);
 
/* Subroutine */ int LAPACKMANGLE(dlasr)(char *side, char *pivot, char *direct, integer *m,
	 integer *n, doublereal *c__, doublereal *s, doublereal *a, integer *
	lda);
 
/* Subroutine */ int LAPACKMANGLE(dlasrt)(char *id, integer *n, doublereal *d__, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(dlassq)(integer *n, doublereal *x, integer *incx, 
	doublereal *scale, doublereal *sumsq);
 
/* Subroutine */ int LAPACKMANGLE(dlasv2)(doublereal *f, doublereal *g, doublereal *h__, 
	doublereal *ssmin, doublereal *ssmax, doublereal *snr, doublereal *
	csr, doublereal *snl, doublereal *csl);
 
/* Subroutine */ int LAPACKMANGLE(dlaswp)(integer *n, doublereal *a, integer *lda, integer 
	*k1, integer *k2, integer *ipiv, integer *incx);
 
/* Subroutine */ int LAPACKMANGLE(dlasy2)(logical *ltranl, logical *ltranr, integer *isgn, 
	integer *n1, integer *n2, doublereal *tl, integer *ldtl, doublereal *
	tr, integer *ldtr, doublereal *b, integer *ldb, doublereal *scale, 
	doublereal *x, integer *ldx, doublereal *xnorm, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dlasyf)(char *uplo, integer *n, integer *nb, integer *kb,
	 doublereal *a, integer *lda, integer *ipiv, doublereal *w, integer *
	ldw, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dlatbs)(char *uplo, char *trans, char *diag, char *
	normin, integer *n, integer *kd, doublereal *ab, integer *ldab, 
	doublereal *x, doublereal *scale, doublereal *cnorm, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dlatdf)(integer *ijob, integer *n, doublereal *z__, 
	integer *ldz, doublereal *rhs, doublereal *rdsum, doublereal *rdscal, 
	integer *ipiv, integer *jpiv);
 
/* Subroutine */ int LAPACKMANGLE(dlatps)(char *uplo, char *trans, char *diag, char *
	normin, integer *n, doublereal *ap, doublereal *x, doublereal *scale, 
	doublereal *cnorm, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dlatrd)(char *uplo, integer *n, integer *nb, doublereal *
	a, integer *lda, doublereal *e, doublereal *tau, doublereal *w, 
	integer *ldw);
 
/* Subroutine */ int LAPACKMANGLE(dlatrs)(char *uplo, char *trans, char *diag, char *
	normin, integer *n, doublereal *a, integer *lda, doublereal *x, 
	doublereal *scale, doublereal *cnorm, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dlatrz)(integer *m, integer *n, integer *l, doublereal *
	a, integer *lda, doublereal *tau, doublereal *work);
 
/* Subroutine */ int LAPACKMANGLE(dlatzm)(char *side, integer *m, integer *n, doublereal *
	v, integer *incv, doublereal *tau, doublereal *c1, doublereal *c2, 
	integer *ldc, doublereal *work);
 
/* Subroutine */ int LAPACKMANGLE(dlauu2)(char *uplo, integer *n, doublereal *a, integer *
	lda, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dlauum)(char *uplo, integer *n, doublereal *a, integer *
	lda, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dopgtr)(char *uplo, integer *n, doublereal *ap, 
	doublereal *tau, doublereal *q, integer *ldq, doublereal *work, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dopmtr)(char *side, char *uplo, char *trans, integer *m, 
	integer *n, doublereal *ap, doublereal *tau, doublereal *c__, integer 
	*ldc, doublereal *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dorg2l)(integer *m, integer *n, integer *k, doublereal *
	a, integer *lda, doublereal *tau, doublereal *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dorg2r)(integer *m, integer *n, integer *k, doublereal *
	a, integer *lda, doublereal *tau, doublereal *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dorgbr)(char *vect, integer *m, integer *n, integer *k, 
	doublereal *a, integer *lda, doublereal *tau, doublereal *work, 
	integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dorghr)(integer *n, integer *ilo, integer *ihi, 
	doublereal *a, integer *lda, doublereal *tau, doublereal *work, 
	integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dorgl2)(integer *m, integer *n, integer *k, doublereal *
	a, integer *lda, doublereal *tau, doublereal *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dorglq)(integer *m, integer *n, integer *k, doublereal *
	a, integer *lda, doublereal *tau, doublereal *work, integer *lwork, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dorgql)(integer *m, integer *n, integer *k, doublereal *
	a, integer *lda, doublereal *tau, doublereal *work, integer *lwork, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dorgqr)(integer *m, integer *n, integer *k, doublereal *
	a, integer *lda, doublereal *tau, doublereal *work, integer *lwork, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dorgr2)(integer *m, integer *n, integer *k, doublereal *
	a, integer *lda, doublereal *tau, doublereal *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dorgrq)(integer *m, integer *n, integer *k, doublereal *
	a, integer *lda, doublereal *tau, doublereal *work, integer *lwork, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dorgtr)(char *uplo, integer *n, doublereal *a, integer *
	lda, doublereal *tau, doublereal *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dorm2l)(char *side, char *trans, integer *m, integer *n, 
	integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *
	c__, integer *ldc, doublereal *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dorm2r)(char *side, char *trans, integer *m, integer *n, 
	integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *
	c__, integer *ldc, doublereal *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dormbr)(char *vect, char *side, char *trans, integer *m, 
	integer *n, integer *k, doublereal *a, integer *lda, doublereal *tau, 
	doublereal *c__, integer *ldc, doublereal *work, integer *lwork, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dormhr)(char *side, char *trans, integer *m, integer *n, 
	integer *ilo, integer *ihi, doublereal *a, integer *lda, doublereal *
	tau, doublereal *c__, integer *ldc, doublereal *work, integer *lwork, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dorml2)(char *side, char *trans, integer *m, integer *n, 
	integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *
	c__, integer *ldc, doublereal *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dormlq)(char *side, char *trans, integer *m, integer *n, 
	integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *
	c__, integer *ldc, doublereal *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dormql)(char *side, char *trans, integer *m, integer *n, 
	integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *
	c__, integer *ldc, doublereal *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dormqr)(char *side, char *trans, integer *m, integer *n, 
	integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *
	c__, integer *ldc, doublereal *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dormr2)(char *side, char *trans, integer *m, integer *n, 
	integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *
	c__, integer *ldc, doublereal *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dormr3)(char *side, char *trans, integer *m, integer *n, 
	integer *k, integer *l, doublereal *a, integer *lda, doublereal *tau, 
	doublereal *c__, integer *ldc, doublereal *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dormrq)(char *side, char *trans, integer *m, integer *n, 
	integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *
	c__, integer *ldc, doublereal *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dormrz)(char *side, char *trans, integer *m, integer *n, 
	integer *k, integer *l, doublereal *a, integer *lda, doublereal *tau, 
	doublereal *c__, integer *ldc, doublereal *work, integer *lwork, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dormtr)(char *side, char *uplo, char *trans, integer *m, 
	integer *n, doublereal *a, integer *lda, doublereal *tau, doublereal *
	c__, integer *ldc, doublereal *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dpbcon)(char *uplo, integer *n, integer *kd, doublereal *
	ab, integer *ldab, doublereal *anorm, doublereal *rcond, doublereal *
	work, integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dpbequ)(char *uplo, integer *n, integer *kd, doublereal *
	ab, integer *ldab, doublereal *s, doublereal *scond, doublereal *amax,
	 integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dpbrfs)(char *uplo, integer *n, integer *kd, integer *
	nrhs, doublereal *ab, integer *ldab, doublereal *afb, integer *ldafb, 
	doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *
	ferr, doublereal *berr, doublereal *work, integer *iwork, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(dpbstf)(char *uplo, integer *n, integer *kd, doublereal *
	ab, integer *ldab, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dpbsv)(char *uplo, integer *n, integer *kd, integer *
	nrhs, doublereal *ab, integer *ldab, doublereal *b, integer *ldb, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dpbsvx)(char *fact, char *uplo, integer *n, integer *kd, 
	integer *nrhs, doublereal *ab, integer *ldab, doublereal *afb, 
	integer *ldafb, char *equed, doublereal *s, doublereal *b, integer *
	ldb, doublereal *x, integer *ldx, doublereal *rcond, doublereal *ferr,
	 doublereal *berr, doublereal *work, integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dpbtf2)(char *uplo, integer *n, integer *kd, doublereal *
	ab, integer *ldab, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dpbtrf)(char *uplo, integer *n, integer *kd, doublereal *
	ab, integer *ldab, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dpbtrs)(char *uplo, integer *n, integer *kd, integer *
	nrhs, doublereal *ab, integer *ldab, doublereal *b, integer *ldb, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dpocon)(char *uplo, integer *n, doublereal *a, integer *
	lda, doublereal *anorm, doublereal *rcond, doublereal *work, integer *
	iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dpoequ)(integer *n, doublereal *a, integer *lda, 
	doublereal *s, doublereal *scond, doublereal *amax, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dporfs)(char *uplo, integer *n, integer *nrhs, 
	doublereal *a, integer *lda, doublereal *af, integer *ldaf, 
	doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *
	ferr, doublereal *berr, doublereal *work, integer *iwork, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(dposv)(char *uplo, integer *n, integer *nrhs, doublereal 
	*a, integer *lda, doublereal *b, integer *ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dposvx)(char *fact, char *uplo, integer *n, integer *
	nrhs, doublereal *a, integer *lda, doublereal *af, integer *ldaf, 
	char *equed, doublereal *s, doublereal *b, integer *ldb, doublereal *
	x, integer *ldx, doublereal *rcond, doublereal *ferr, doublereal *
	berr, doublereal *work, integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dpotf2)(char *uplo, integer *n, doublereal *a, integer *
	lda, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dpotrf)(char *uplo, integer *n, doublereal *a, integer *
	lda, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dpotri)(char *uplo, integer *n, doublereal *a, integer *
	lda, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dpotrs)(char *uplo, integer *n, integer *nrhs, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(dppcon)(char *uplo, integer *n, doublereal *ap, 
	doublereal *anorm, doublereal *rcond, doublereal *work, integer *
	iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dppequ)(char *uplo, integer *n, doublereal *ap, 
	doublereal *s, doublereal *scond, doublereal *amax, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dpprfs)(char *uplo, integer *n, integer *nrhs, 
	doublereal *ap, doublereal *afp, doublereal *b, integer *ldb, 
	doublereal *x, integer *ldx, doublereal *ferr, doublereal *berr, 
	doublereal *work, integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dppsv)(char *uplo, integer *n, integer *nrhs, doublereal 
	*ap, doublereal *b, integer *ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dppsvx)(char *fact, char *uplo, integer *n, integer *
	nrhs, doublereal *ap, doublereal *afp, char *equed, doublereal *s, 
	doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *
	rcond, doublereal *ferr, doublereal *berr, doublereal *work, integer *
	iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dpptrf)(char *uplo, integer *n, doublereal *ap, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(dpptri)(char *uplo, integer *n, doublereal *ap, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(dpptrs)(char *uplo, integer *n, integer *nrhs, 
	doublereal *ap, doublereal *b, integer *ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dptcon)(integer *n, doublereal *d__, doublereal *e, 
	doublereal *anorm, doublereal *rcond, doublereal *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dpteqr)(char *compz, integer *n, doublereal *d__, 
	doublereal *e, doublereal *z__, integer *ldz, doublereal *work, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dptrfs)(integer *n, integer *nrhs, doublereal *d__, 
	doublereal *e, doublereal *df, doublereal *ef, doublereal *b, integer 
	*ldb, doublereal *x, integer *ldx, doublereal *ferr, doublereal *berr,
	 doublereal *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dptsv)(integer *n, integer *nrhs, doublereal *d__, 
	doublereal *e, doublereal *b, integer *ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dptsvx)(char *fact, integer *n, integer *nrhs, 
	doublereal *d__, doublereal *e, doublereal *df, doublereal *ef, 
	doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *
	rcond, doublereal *ferr, doublereal *berr, doublereal *work, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(dpttrf)(integer *n, doublereal *d__, doublereal *e, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dpttrs)(integer *n, integer *nrhs, doublereal *d__, 
	doublereal *e, doublereal *b, integer *ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dptts2)(integer *n, integer *nrhs, doublereal *d__, 
	doublereal *e, doublereal *b, integer *ldb);
 
/* Subroutine */ int LAPACKMANGLE(drscl)(integer *n, doublereal *sa, doublereal *sx, 
	integer *incx);
 
/* Subroutine */ int LAPACKMANGLE(dsbev)(char *jobz, char *uplo, integer *n, integer *kd, 
	doublereal *ab, integer *ldab, doublereal *w, doublereal *z__, 
	integer *ldz, doublereal *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dsbevd)(char *jobz, char *uplo, integer *n, integer *kd, 
	doublereal *ab, integer *ldab, doublereal *w, doublereal *z__, 
	integer *ldz, doublereal *work, integer *lwork, integer *iwork, 
	integer *liwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dsbevx)(char *jobz, char *range, char *uplo, integer *n, 
	integer *kd, doublereal *ab, integer *ldab, doublereal *q, integer *
	ldq, doublereal *vl, doublereal *vu, integer *il, integer *iu, 
	doublereal *abstol, integer *m, doublereal *w, doublereal *z__, 
	integer *ldz, doublereal *work, integer *iwork, integer *ifail, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dsbgst)(char *vect, char *uplo, integer *n, integer *ka, 
	integer *kb, doublereal *ab, integer *ldab, doublereal *bb, integer *
	ldbb, doublereal *x, integer *ldx, doublereal *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dsbgv)(char *jobz, char *uplo, integer *n, integer *ka, 
	integer *kb, doublereal *ab, integer *ldab, doublereal *bb, integer *
	ldbb, doublereal *w, doublereal *z__, integer *ldz, doublereal *work, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dsbgvd)(char *jobz, char *uplo, integer *n, integer *ka, 
	integer *kb, doublereal *ab, integer *ldab, doublereal *bb, integer *
	ldbb, doublereal *w, doublereal *z__, integer *ldz, doublereal *work, 
	integer *lwork, integer *iwork, integer *liwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dsbgvx)(char *jobz, char *range, char *uplo, integer *n, 
	integer *ka, integer *kb, doublereal *ab, integer *ldab, doublereal *
	bb, integer *ldbb, doublereal *q, integer *ldq, doublereal *vl, 
	doublereal *vu, integer *il, integer *iu, doublereal *abstol, integer 
	*m, doublereal *w, doublereal *z__, integer *ldz, doublereal *work, 
	integer *iwork, integer *ifail, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dsbtrd)(char *vect, char *uplo, integer *n, integer *kd, 
	doublereal *ab, integer *ldab, doublereal *d__, doublereal *e, 
	doublereal *q, integer *ldq, doublereal *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dspcon)(char *uplo, integer *n, doublereal *ap, integer *
	ipiv, doublereal *anorm, doublereal *rcond, doublereal *work, integer 
	*iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dspev)(char *jobz, char *uplo, integer *n, doublereal *
	ap, doublereal *w, doublereal *z__, integer *ldz, doublereal *work, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dspevd)(char *jobz, char *uplo, integer *n, doublereal *
	ap, doublereal *w, doublereal *z__, integer *ldz, doublereal *work, 
	integer *lwork, integer *iwork, integer *liwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dspevx)(char *jobz, char *range, char *uplo, integer *n, 
	doublereal *ap, doublereal *vl, doublereal *vu, integer *il, integer *
	iu, doublereal *abstol, integer *m, doublereal *w, doublereal *z__, 
	integer *ldz, doublereal *work, integer *iwork, integer *ifail, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dspgst)(integer *itype, char *uplo, integer *n, 
	doublereal *ap, doublereal *bp, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dspgv)(integer *itype, char *jobz, char *uplo, integer *
	n, doublereal *ap, doublereal *bp, doublereal *w, doublereal *z__, 
	integer *ldz, doublereal *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dspgvd)(integer *itype, char *jobz, char *uplo, integer *
	n, doublereal *ap, doublereal *bp, doublereal *w, doublereal *z__, 
	integer *ldz, doublereal *work, integer *lwork, integer *iwork, 
	integer *liwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dspgvx)(integer *itype, char *jobz, char *range, char *
	uplo, integer *n, doublereal *ap, doublereal *bp, doublereal *vl, 
	doublereal *vu, integer *il, integer *iu, doublereal *abstol, integer 
	*m, doublereal *w, doublereal *z__, integer *ldz, doublereal *work, 
	integer *iwork, integer *ifail, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dsprfs)(char *uplo, integer *n, integer *nrhs, 
	doublereal *ap, doublereal *afp, integer *ipiv, doublereal *b, 
	integer *ldb, doublereal *x, integer *ldx, doublereal *ferr, 
	doublereal *berr, doublereal *work, integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dspsv)(char *uplo, integer *n, integer *nrhs, doublereal 
	*ap, integer *ipiv, doublereal *b, integer *ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dspsvx)(char *fact, char *uplo, integer *n, integer *
	nrhs, doublereal *ap, doublereal *afp, integer *ipiv, doublereal *b, 
	integer *ldb, doublereal *x, integer *ldx, doublereal *rcond, 
	doublereal *ferr, doublereal *berr, doublereal *work, integer *iwork, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dsptrd)(char *uplo, integer *n, doublereal *ap, 
	doublereal *d__, doublereal *e, doublereal *tau, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dsptrf)(char *uplo, integer *n, doublereal *ap, integer *
	ipiv, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dsptri)(char *uplo, integer *n, doublereal *ap, integer *
	ipiv, doublereal *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dsptrs)(char *uplo, integer *n, integer *nrhs, 
	doublereal *ap, integer *ipiv, doublereal *b, integer *ldb, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(dstebz)(char *range, char *order, integer *n, doublereal 
	*vl, doublereal *vu, integer *il, integer *iu, doublereal *abstol, 
	doublereal *d__, doublereal *e, integer *m, integer *nsplit, 
	doublereal *w, integer *iblock, integer *isplit, doublereal *work, 
	integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dstedc)(char *compz, integer *n, doublereal *d__, 
	doublereal *e, doublereal *z__, integer *ldz, doublereal *work, 
	integer *lwork, integer *iwork, integer *liwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dstegr)(char *jobz, char *range, integer *n, doublereal *
	d__, doublereal *e, doublereal *vl, doublereal *vu, integer *il, 
	integer *iu, doublereal *abstol, integer *m, doublereal *w, 
	doublereal *z__, integer *ldz, integer *isuppz, doublereal *work, 
	integer *lwork, integer *iwork, integer *liwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dstein)(integer *n, doublereal *d__, doublereal *e, 
	integer *m, doublereal *w, integer *iblock, integer *isplit, 
	doublereal *z__, integer *ldz, doublereal *work, integer *iwork, 
	integer *ifail, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dsteqr)(char *compz, integer *n, doublereal *d__, 
	doublereal *e, doublereal *z__, integer *ldz, doublereal *work, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dsterf)(integer *n, doublereal *d__, doublereal *e, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dstev)(char *jobz, integer *n, doublereal *d__, 
	doublereal *e, doublereal *z__, integer *ldz, doublereal *work, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dstevd)(char *jobz, integer *n, doublereal *d__, 
	doublereal *e, doublereal *z__, integer *ldz, doublereal *work, 
	integer *lwork, integer *iwork, integer *liwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dstevr)(char *jobz, char *range, integer *n, doublereal *
	d__, doublereal *e, doublereal *vl, doublereal *vu, integer *il, 
	integer *iu, doublereal *abstol, integer *m, doublereal *w, 
	doublereal *z__, integer *ldz, integer *isuppz, doublereal *work, 
	integer *lwork, integer *iwork, integer *liwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dstevx)(char *jobz, char *range, integer *n, doublereal *
	d__, doublereal *e, doublereal *vl, doublereal *vu, integer *il, 
	integer *iu, doublereal *abstol, integer *m, doublereal *w, 
	doublereal *z__, integer *ldz, doublereal *work, integer *iwork, 
	integer *ifail, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dsycon)(char *uplo, integer *n, doublereal *a, integer *
	lda, integer *ipiv, doublereal *anorm, doublereal *rcond, doublereal *
	work, integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dsyev)(char *jobz, char *uplo, integer *n, doublereal *a,
	 integer *lda, doublereal *w, doublereal *work, integer *lwork, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dsyevd)(char *jobz, char *uplo, integer *n, doublereal *
	a, integer *lda, doublereal *w, doublereal *work, integer *lwork, 
	integer *iwork, integer *liwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dsyevr)(char *jobz, char *range, char *uplo, integer *n, 
	doublereal *a, integer *lda, doublereal *vl, doublereal *vu, integer *
	il, integer *iu, doublereal *abstol, integer *m, doublereal *w, 
	doublereal *z__, integer *ldz, integer *isuppz, doublereal *work, 
	integer *lwork, integer *iwork, integer *liwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dsyevx)(char *jobz, char *range, char *uplo, integer *n, 
	doublereal *a, integer *lda, doublereal *vl, doublereal *vu, integer *
	il, integer *iu, doublereal *abstol, integer *m, doublereal *w, 
	doublereal *z__, integer *ldz, doublereal *work, integer *lwork, 
	integer *iwork, integer *ifail, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dsygs2)(integer *itype, char *uplo, integer *n, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(dsygst)(integer *itype, char *uplo, integer *n, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(dsygv)(integer *itype, char *jobz, char *uplo, integer *
	n, doublereal *a, integer *lda, doublereal *b, integer *ldb, 
	doublereal *w, doublereal *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dsygvd)(integer *itype, char *jobz, char *uplo, integer *
	n, doublereal *a, integer *lda, doublereal *b, integer *ldb, 
	doublereal *w, doublereal *work, integer *lwork, integer *iwork, 
	integer *liwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dsygvx)(integer *itype, char *jobz, char *range, char *
	uplo, integer *n, doublereal *a, integer *lda, doublereal *b, integer 
	*ldb, doublereal *vl, doublereal *vu, integer *il, integer *iu, 
	doublereal *abstol, integer *m, doublereal *w, doublereal *z__, 
	integer *ldz, doublereal *work, integer *lwork, integer *iwork, 
	integer *ifail, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dsyrfs)(char *uplo, integer *n, integer *nrhs, 
	doublereal *a, integer *lda, doublereal *af, integer *ldaf, integer *
	ipiv, doublereal *b, integer *ldb, doublereal *x, integer *ldx, 
	doublereal *ferr, doublereal *berr, doublereal *work, integer *iwork, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dsysv)(char *uplo, integer *n, integer *nrhs, doublereal 
	*a, integer *lda, integer *ipiv, doublereal *b, integer *ldb, 
	doublereal *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dsysvx)(char *fact, char *uplo, integer *n, integer *
	nrhs, doublereal *a, integer *lda, doublereal *af, integer *ldaf, 
	integer *ipiv, doublereal *b, integer *ldb, doublereal *x, integer *
	ldx, doublereal *rcond, doublereal *ferr, doublereal *berr, 
	doublereal *work, integer *lwork, integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dsytd2)(char *uplo, integer *n, doublereal *a, integer *
	lda, doublereal *d__, doublereal *e, doublereal *tau, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dsytf2)(char *uplo, integer *n, doublereal *a, integer *
	lda, integer *ipiv, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dsytrd)(char *uplo, integer *n, doublereal *a, integer *
	lda, doublereal *d__, doublereal *e, doublereal *tau, doublereal *
	work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dsytrf)(char *uplo, integer *n, doublereal *a, integer *
	lda, integer *ipiv, doublereal *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dsytri)(char *uplo, integer *n, doublereal *a, integer *
	lda, integer *ipiv, doublereal *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dsytrs)(char *uplo, integer *n, integer *nrhs, 
	doublereal *a, integer *lda, integer *ipiv, doublereal *b, integer *
	ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dtbcon)(char *norm, char *uplo, char *diag, integer *n, 
	integer *kd, doublereal *ab, integer *ldab, doublereal *rcond, 
	doublereal *work, integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dtbrfs)(char *uplo, char *trans, char *diag, integer *n, 
	integer *kd, integer *nrhs, doublereal *ab, integer *ldab, doublereal 
	*b, integer *ldb, doublereal *x, integer *ldx, doublereal *ferr, 
	doublereal *berr, doublereal *work, integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dtbtrs)(char *uplo, char *trans, char *diag, integer *n, 
	integer *kd, integer *nrhs, doublereal *ab, integer *ldab, doublereal 
	*b, integer *ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dtgevc)(char *side, char *howmny, logical *select, 
	integer *n, doublereal *a, integer *lda, doublereal *b, integer *ldb, 
	doublereal *vl, integer *ldvl, doublereal *vr, integer *ldvr, integer 
	*mm, integer *m, doublereal *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dtgex2)(logical *wantq, logical *wantz, integer *n, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	q, integer *ldq, doublereal *z__, integer *ldz, integer *j1, integer *
	n1, integer *n2, doublereal *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dtgexc)(logical *wantq, logical *wantz, integer *n, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	q, integer *ldq, doublereal *z__, integer *ldz, integer *ifst, 
	integer *ilst, doublereal *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dtgsen)(integer *ijob, logical *wantq, logical *wantz, 
	logical *select, integer *n, doublereal *a, integer *lda, doublereal *
	b, integer *ldb, doublereal *alphar, doublereal *alphai, doublereal *
	beta, doublereal *q, integer *ldq, doublereal *z__, integer *ldz, 
	integer *m, doublereal *pl, doublereal *pr, doublereal *dif, 
	doublereal *work, integer *lwork, integer *iwork, integer *liwork, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dtgsja)(char *jobu, char *jobv, char *jobq, integer *m, 
	integer *p, integer *n, integer *k, integer *l, doublereal *a, 
	integer *lda, doublereal *b, integer *ldb, doublereal *tola, 
	doublereal *tolb, doublereal *alpha, doublereal *beta, doublereal *u, 
	integer *ldu, doublereal *v, integer *ldv, doublereal *q, integer *
	ldq, doublereal *work, integer *ncycle, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dtgsna)(char *job, char *howmny, logical *select, 
	integer *n, doublereal *a, integer *lda, doublereal *b, integer *ldb, 
	doublereal *vl, integer *ldvl, doublereal *vr, integer *ldvr, 
	doublereal *s, doublereal *dif, integer *mm, integer *m, doublereal *
	work, integer *lwork, integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dtgsy2)(char *trans, integer *ijob, integer *m, integer *
	n, doublereal *a, integer *lda, doublereal *b, integer *ldb, 
	doublereal *c__, integer *ldc, doublereal *d__, integer *ldd, 
	doublereal *e, integer *lde, doublereal *f, integer *ldf, doublereal *
	scale, doublereal *rdsum, doublereal *rdscal, integer *iwork, integer 
	*pq, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dtgsyl)(char *trans, integer *ijob, integer *m, integer *
	n, doublereal *a, integer *lda, doublereal *b, integer *ldb, 
	doublereal *c__, integer *ldc, doublereal *d__, integer *ldd, 
	doublereal *e, integer *lde, doublereal *f, integer *ldf, doublereal *
	scale, doublereal *dif, doublereal *work, integer *lwork, integer *
	iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dtpcon)(char *norm, char *uplo, char *diag, integer *n, 
	doublereal *ap, doublereal *rcond, doublereal *work, integer *iwork, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dtprfs)(char *uplo, char *trans, char *diag, integer *n, 
	integer *nrhs, doublereal *ap, doublereal *b, integer *ldb, 
	doublereal *x, integer *ldx, doublereal *ferr, doublereal *berr, 
	doublereal *work, integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dtptri)(char *uplo, char *diag, integer *n, doublereal *
	ap, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dtptrs)(char *uplo, char *trans, char *diag, integer *n, 
	integer *nrhs, doublereal *ap, doublereal *b, integer *ldb, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(dtrcon)(char *norm, char *uplo, char *diag, integer *n, 
	doublereal *a, integer *lda, doublereal *rcond, doublereal *work, 
	integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dtrevc)(char *side, char *howmny, logical *select, 
	integer *n, doublereal *t, integer *ldt, doublereal *vl, integer *
	ldvl, doublereal *vr, integer *ldvr, integer *mm, integer *m, 
	doublereal *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dtrexc)(char *compq, integer *n, doublereal *t, integer *
	ldt, doublereal *q, integer *ldq, integer *ifst, integer *ilst, 
	doublereal *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dtrrfs)(char *uplo, char *trans, char *diag, integer *n, 
	integer *nrhs, doublereal *a, integer *lda, doublereal *b, integer *
	ldb, doublereal *x, integer *ldx, doublereal *ferr, doublereal *berr, 
	doublereal *work, integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dtrsen)(char *job, char *compq, logical *select, integer 
	*n, doublereal *t, integer *ldt, doublereal *q, integer *ldq, 
	doublereal *wr, doublereal *wi, integer *m, doublereal *s, doublereal 
	*sep, doublereal *work, integer *lwork, integer *iwork, integer *
	liwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dtrsna)(char *job, char *howmny, logical *select, 
	integer *n, doublereal *t, integer *ldt, doublereal *vl, integer *
	ldvl, doublereal *vr, integer *ldvr, doublereal *s, doublereal *sep, 
	integer *mm, integer *m, doublereal *work, integer *ldwork, integer *
	iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dtrsyl)(char *trana, char *tranb, integer *isgn, integer 
	*m, integer *n, doublereal *a, integer *lda, doublereal *b, integer *
	ldb, doublereal *c__, integer *ldc, doublereal *scale, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dtrti2)(char *uplo, char *diag, integer *n, doublereal *
	a, integer *lda, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dtrtri)(char *uplo, char *diag, integer *n, doublereal *
	a, integer *lda, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dtrtrs)(char *uplo, char *trans, char *diag, integer *n, 
	integer *nrhs, doublereal *a, integer *lda, doublereal *b, integer *
	ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dtzrqf)(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *tau, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(dtzrzf)(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *tau, doublereal *work, integer *lwork, integer *info);
 
integer icmax1_(integer *n, complex *cx, integer *incx);
 
integer ieeeck_(integer *ispec, real *zero, real *one);
 
integer ilaenv_(integer *ispec, char *name__, char *opts, integer *n1, 
	integer *n2, integer *n3, integer *n4, ftnlen name_len, ftnlen 
	opts_len);
 
integer izmax1_(integer *n, doublecomplex *cx, integer *incx);
 
/* Subroutine */ int LAPACKMANGLE(sbdsdc)(char *uplo, char *compq, integer *n, real *d__, 
	real *e, real *u, integer *ldu, real *vt, integer *ldvt, real *q, 
	integer *iq, real *work, integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sbdsqr)(char *uplo, integer *n, integer *ncvt, integer *
	nru, integer *ncc, real *d__, real *e, real *vt, integer *ldvt, real *
	u, integer *ldu, real *c__, integer *ldc, real *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sdisna)(char *job, integer *m, integer *n, real *d__, 
	real *sep, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sgbbrd)(char *vect, integer *m, integer *n, integer *ncc,
	 integer *kl, integer *ku, real *ab, integer *ldab, real *d__, real *
	e, real *q, integer *ldq, real *pt, integer *ldpt, real *c__, integer 
	*ldc, real *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sgbcon)(char *norm, integer *n, integer *kl, integer *ku,
	 real *ab, integer *ldab, integer *ipiv, real *anorm, real *rcond, 
	real *work, integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sgbequ)(integer *m, integer *n, integer *kl, integer *ku,
	 real *ab, integer *ldab, real *r__, real *c__, real *rowcnd, real *
	colcnd, real *amax, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sgbrfs)(char *trans, integer *n, integer *kl, integer *
	ku, integer *nrhs, real *ab, integer *ldab, real *afb, integer *ldafb,
	 integer *ipiv, real *b, integer *ldb, real *x, integer *ldx, real *
	ferr, real *berr, real *work, integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sgbsv)(integer *n, integer *kl, integer *ku, integer *
	nrhs, real *ab, integer *ldab, integer *ipiv, real *b, integer *ldb, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sgbsvx)(char *fact, char *trans, integer *n, integer *kl,
	 integer *ku, integer *nrhs, real *ab, integer *ldab, real *afb, 
	integer *ldafb, integer *ipiv, char *equed, real *r__, real *c__, 
	real *b, integer *ldb, real *x, integer *ldx, real *rcond, real *ferr,
	 real *berr, real *work, integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sgbtf2)(integer *m, integer *n, integer *kl, integer *ku,
	 real *ab, integer *ldab, integer *ipiv, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sgbtrf)(integer *m, integer *n, integer *kl, integer *ku,
	 real *ab, integer *ldab, integer *ipiv, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sgbtrs)(char *trans, integer *n, integer *kl, integer *
	ku, integer *nrhs, real *ab, integer *ldab, integer *ipiv, real *b, 
	integer *ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sgebak)(char *job, char *side, integer *n, integer *ilo, 
	integer *ihi, real *scale, integer *m, real *v, integer *ldv, integer 
	*info);
 
/* Subroutine */ int LAPACKMANGLE(sgebal)(char *job, integer *n, real *a, integer *lda, 
	integer *ilo, integer *ihi, real *scale, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sgebd2)(integer *m, integer *n, real *a, integer *lda, 
	real *d__, real *e, real *tauq, real *taup, real *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sgebrd)(integer *m, integer *n, real *a, integer *lda, 
	real *d__, real *e, real *tauq, real *taup, real *work, integer *
	lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sgecon)(char *norm, integer *n, real *a, integer *lda, 
	real *anorm, real *rcond, real *work, integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sgeequ)(integer *m, integer *n, real *a, integer *lda, 
	real *r__, real *c__, real *rowcnd, real *colcnd, real *amax, integer 
	*info);
 
/* Subroutine */ int LAPACKMANGLE(sgees)(char *jobvs, char *sort, L_fp select, integer *n, 
	real *a, integer *lda, integer *sdim, real *wr, real *wi, real *vs, 
	integer *ldvs, real *work, integer *lwork, logical *bwork, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(sgeesx)(char *jobvs, char *sort, L_fp select, char *
	sense, integer *n, real *a, integer *lda, integer *sdim, real *wr, 
	real *wi, real *vs, integer *ldvs, real *rconde, real *rcondv, real *
	work, integer *lwork, integer *iwork, integer *liwork, logical *bwork,
	 integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sgeev)(char *jobvl, char *jobvr, integer *n, real *a, 
	integer *lda, real *wr, real *wi, real *vl, integer *ldvl, real *vr, 
	integer *ldvr, real *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sgeevx)(char *balanc, char *jobvl, char *jobvr, char *
	sense, integer *n, real *a, integer *lda, real *wr, real *wi, real *
	vl, integer *ldvl, real *vr, integer *ldvr, integer *ilo, integer *
	ihi, real *scale, real *abnrm, real *rconde, real *rcondv, real *work,
	 integer *lwork, integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sgegs)(char *jobvsl, char *jobvsr, integer *n, real *a, 
	integer *lda, real *b, integer *ldb, real *alphar, real *alphai, real 
	*beta, real *vsl, integer *ldvsl, real *vsr, integer *ldvsr, real *
	work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sgegv)(char *jobvl, char *jobvr, integer *n, real *a, 
	integer *lda, real *b, integer *ldb, real *alphar, real *alphai, real 
	*beta, real *vl, integer *ldvl, real *vr, integer *ldvr, real *work, 
	integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sgehd2)(integer *n, integer *ilo, integer *ihi, real *a, 
	integer *lda, real *tau, real *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sgehrd)(integer *n, integer *ilo, integer *ihi, real *a, 
	integer *lda, real *tau, real *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sgelq2)(integer *m, integer *n, real *a, integer *lda, 
	real *tau, real *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sgelqf)(integer *m, integer *n, real *a, integer *lda, 
	real *tau, real *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sgels)(char *trans, integer *m, integer *n, integer *
	nrhs, real *a, integer *lda, real *b, integer *ldb, real *work, 
	integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sgelsd)(integer *m, integer *n, integer *nrhs, real *a, 
	integer *lda, real *b, integer *ldb, real *s, real *rcond, integer *
	rank, real *work, integer *lwork, integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sgelss)(integer *m, integer *n, integer *nrhs, real *a, 
	integer *lda, real *b, integer *ldb, real *s, real *rcond, integer *
	rank, real *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sgelsx)(integer *m, integer *n, integer *nrhs, real *a, 
	integer *lda, real *b, integer *ldb, integer *jpvt, real *rcond, 
	integer *rank, real *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sgelsy)(integer *m, integer *n, integer *nrhs, real *a, 
	integer *lda, real *b, integer *ldb, integer *jpvt, real *rcond, 
	integer *rank, real *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sgeql2)(integer *m, integer *n, real *a, integer *lda, 
	real *tau, real *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sgeqlf)(integer *m, integer *n, real *a, integer *lda, 
	real *tau, real *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sgeqp3)(integer *m, integer *n, real *a, integer *lda, 
	integer *jpvt, real *tau, real *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sgeqpf)(integer *m, integer *n, real *a, integer *lda, 
	integer *jpvt, real *tau, real *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sgeqr2)(integer *m, integer *n, real *a, integer *lda, 
	real *tau, real *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sgeqrf)(integer *m, integer *n, real *a, integer *lda, 
	real *tau, real *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sgerfs)(char *trans, integer *n, integer *nrhs, real *a, 
	integer *lda, real *af, integer *ldaf, integer *ipiv, real *b, 
	integer *ldb, real *x, integer *ldx, real *ferr, real *berr, real *
	work, integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sgerq2)(integer *m, integer *n, real *a, integer *lda, 
	real *tau, real *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sgerqf)(integer *m, integer *n, real *a, integer *lda, 
	real *tau, real *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sgesc2)(integer *n, real *a, integer *lda, real *rhs, 
	integer *ipiv, integer *jpiv, real *scale);
 
/* Subroutine */ int LAPACKMANGLE(sgesdd)(char *jobz, integer *m, integer *n, real *a, 
	integer *lda, real *s, real *u, integer *ldu, real *vt, integer *ldvt,
	 real *work, integer *lwork, integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sgesv)(integer *n, integer *nrhs, real *a, integer *lda, 
	integer *ipiv, real *b, integer *ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sgesvd)(char *jobu, char *jobvt, integer *m, integer *n, 
	real *a, integer *lda, real *s, real *u, integer *ldu, real *vt, 
	integer *ldvt, real *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sgesvx)(char *fact, char *trans, integer *n, integer *
	nrhs, real *a, integer *lda, real *af, integer *ldaf, integer *ipiv, 
	char *equed, real *r__, real *c__, real *b, integer *ldb, real *x, 
	integer *ldx, real *rcond, real *ferr, real *berr, real *work, 
	integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sgetc2)(integer *n, real *a, integer *lda, integer *ipiv,
	 integer *jpiv, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sgetf2)(integer *m, integer *n, real *a, integer *lda, 
	integer *ipiv, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sgetrf)(integer *m, integer *n, real *a, integer *lda, 
	integer *ipiv, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sgetri)(integer *n, real *a, integer *lda, integer *ipiv,
	 real *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sgetrs)(char *trans, integer *n, integer *nrhs, real *a, 
	integer *lda, integer *ipiv, real *b, integer *ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sggbak)(char *job, char *side, integer *n, integer *ilo, 
	integer *ihi, real *lscale, real *rscale, integer *m, real *v, 
	integer *ldv, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sggbal)(char *job, integer *n, real *a, integer *lda, 
	real *b, integer *ldb, integer *ilo, integer *ihi, real *lscale, real 
	*rscale, real *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sgges)(char *jobvsl, char *jobvsr, char *sort, L_fp 
	selctg, integer *n, real *a, integer *lda, real *b, integer *ldb, 
	integer *sdim, real *alphar, real *alphai, real *beta, real *vsl, 
	integer *ldvsl, real *vsr, integer *ldvsr, real *work, integer *lwork,
	 logical *bwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sggesx)(char *jobvsl, char *jobvsr, char *sort, L_fp 
	selctg, char *sense, integer *n, real *a, integer *lda, real *b, 
	integer *ldb, integer *sdim, real *alphar, real *alphai, real *beta, 
	real *vsl, integer *ldvsl, real *vsr, integer *ldvsr, real *rconde, 
	real *rcondv, real *work, integer *lwork, integer *iwork, integer *
	liwork, logical *bwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sggev)(char *jobvl, char *jobvr, integer *n, real *a, 
	integer *lda, real *b, integer *ldb, real *alphar, real *alphai, real 
	*beta, real *vl, integer *ldvl, real *vr, integer *ldvr, real *work, 
	integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sggevx)(char *balanc, char *jobvl, char *jobvr, char *
	sense, integer *n, real *a, integer *lda, real *b, integer *ldb, real 
	*alphar, real *alphai, real *beta, real *vl, integer *ldvl, real *vr, 
	integer *ldvr, integer *ilo, integer *ihi, real *lscale, real *rscale,
	 real *abnrm, real *bbnrm, real *rconde, real *rcondv, real *work, 
	integer *lwork, integer *iwork, logical *bwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sggglm)(integer *n, integer *m, integer *p, real *a, 
	integer *lda, real *b, integer *ldb, real *d__, real *x, real *y, 
	real *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sgghrd)(char *compq, char *compz, integer *n, integer *
	ilo, integer *ihi, real *a, integer *lda, real *b, integer *ldb, real 
	*q, integer *ldq, real *z__, integer *ldz, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sgglse)(integer *m, integer *n, integer *p, real *a, 
	integer *lda, real *b, integer *ldb, real *c__, real *d__, real *x, 
	real *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sggqrf)(integer *n, integer *m, integer *p, real *a, 
	integer *lda, real *taua, real *b, integer *ldb, real *taub, real *
	work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sggrqf)(integer *m, integer *p, integer *n, real *a, 
	integer *lda, real *taua, real *b, integer *ldb, real *taub, real *
	work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sggsvd)(char *jobu, char *jobv, char *jobq, integer *m, 
	integer *n, integer *p, integer *k, integer *l, real *a, integer *lda,
	 real *b, integer *ldb, real *alpha, real *beta, real *u, integer *
	ldu, real *v, integer *ldv, real *q, integer *ldq, real *work, 
	integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sggsvp)(char *jobu, char *jobv, char *jobq, integer *m, 
	integer *p, integer *n, real *a, integer *lda, real *b, integer *ldb, 
	real *tola, real *tolb, integer *k, integer *l, real *u, integer *ldu,
	 real *v, integer *ldv, real *q, integer *ldq, integer *iwork, real *
	tau, real *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sgtcon)(char *norm, integer *n, real *dl, real *d__, 
	real *du, real *du2, integer *ipiv, real *anorm, real *rcond, real *
	work, integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sgtrfs)(char *trans, integer *n, integer *nrhs, real *dl,
	 real *d__, real *du, real *dlf, real *df, real *duf, real *du2, 
	integer *ipiv, real *b, integer *ldb, real *x, integer *ldx, real *
	ferr, real *berr, real *work, integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sgtsv)(integer *n, integer *nrhs, real *dl, real *d__, 
	real *du, real *b, integer *ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sgtsvx)(char *fact, char *trans, integer *n, integer *
	nrhs, real *dl, real *d__, real *du, real *dlf, real *df, real *duf, 
	real *du2, integer *ipiv, real *b, integer *ldb, real *x, integer *
	ldx, real *rcond, real *ferr, real *berr, real *work, integer *iwork, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sgttrf)(integer *n, real *dl, real *d__, real *du, real *
	du2, integer *ipiv, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sgttrs)(char *trans, integer *n, integer *nrhs, real *dl,
	 real *d__, real *du, real *du2, integer *ipiv, real *b, integer *ldb,
	 integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sgtts2)(integer *itrans, integer *n, integer *nrhs, real 
	*dl, real *d__, real *du, real *du2, integer *ipiv, real *b, integer *
	ldb);
 
/* Subroutine */ int LAPACKMANGLE(shgeqz)(char *job, char *compq, char *compz, integer *n, 
	integer *ilo, integer *ihi, real *a, integer *lda, real *b, integer *
	ldb, real *alphar, real *alphai, real *beta, real *q, integer *ldq, 
	real *z__, integer *ldz, real *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(shsein)(char *side, char *eigsrc, char *initv, logical *
	select, integer *n, real *h__, integer *ldh, real *wr, real *wi, real 
	*vl, integer *ldvl, real *vr, integer *ldvr, integer *mm, integer *m, 
	real *work, integer *ifaill, integer *ifailr, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(shseqr)(char *job, char *compz, integer *n, integer *ilo,
	 integer *ihi, real *h__, integer *ldh, real *wr, real *wi, real *z__,
	 integer *ldz, real *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(slabad)(real *small, real *large);
 
/* Subroutine */ int LAPACKMANGLE(slabrd)(integer *m, integer *n, integer *nb, real *a, 
	integer *lda, real *d__, real *e, real *tauq, real *taup, real *x, 
	integer *ldx, real *y, integer *ldy);
 
/* Subroutine */ int LAPACKMANGLE(slacon)(integer *n, real *v, real *x, integer *isgn, 
	real *est, integer *kase);
 
/* Subroutine */ int LAPACKMANGLE(slacpy)(char *uplo, integer *m, integer *n, real *a, 
	integer *lda, real *b, integer *ldb);
 
/* Subroutine */ int LAPACKMANGLE(sladiv)(real *a, real *b, real *c__, real *d__, real *p, 
	real *q);
 
/* Subroutine */ int LAPACKMANGLE(slae2)(real *a, real *b, real *c__, real *rt1, real *rt2);
 
/* Subroutine */ int LAPACKMANGLE(slaebz)(integer *ijob, integer *nitmax, integer *n, 
	integer *mmax, integer *minp, integer *nbmin, real *abstol, real *
	reltol, real *pivmin, real *d__, real *e, real *e2, integer *nval, 
	real *ab, real *c__, integer *mout, integer *nab, real *work, integer 
	*iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(slaed0)(integer *icompq, integer *qsiz, integer *n, real 
	*d__, real *e, real *q, integer *ldq, real *qstore, integer *ldqs, 
	real *work, integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(slaed1)(integer *n, real *d__, real *q, integer *ldq, 
	integer *indxq, real *rho, integer *cutpnt, real *work, integer *
	iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(slaed2)(integer *k, integer *n, integer *n1, real *d__, 
	real *q, integer *ldq, integer *indxq, real *rho, real *z__, real *
	dlamda, real *w, real *q2, integer *indx, integer *indxc, integer *
	indxp, integer *coltyp, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(slaed3)(integer *k, integer *n, integer *n1, real *d__, 
	real *q, integer *ldq, real *rho, real *dlamda, real *q2, integer *
	indx, integer *ctot, real *w, real *s, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(slaed4)(integer *n, integer *i__, real *d__, real *z__, 
	real *delta, real *rho, real *dlam, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(slaed5)(integer *i__, real *d__, real *z__, real *delta, 
	real *rho, real *dlam);
 
/* Subroutine */ int LAPACKMANGLE(slaed6)(integer *kniter, logical *orgati, real *rho, 
	real *d__, real *z__, real *finit, real *tau, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(slaed7)(integer *icompq, integer *n, integer *qsiz, 
	integer *tlvls, integer *curlvl, integer *curpbm, real *d__, real *q, 
	integer *ldq, integer *indxq, real *rho, integer *cutpnt, real *
	qstore, integer *qptr, integer *prmptr, integer *perm, integer *
	givptr, integer *givcol, real *givnum, real *work, integer *iwork, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(slaed8)(integer *icompq, integer *k, integer *n, integer 
	*qsiz, real *d__, real *q, integer *ldq, integer *indxq, real *rho, 
	integer *cutpnt, real *z__, real *dlamda, real *q2, integer *ldq2, 
	real *w, integer *perm, integer *givptr, integer *givcol, real *
	givnum, integer *indxp, integer *indx, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(slaed9)(integer *k, integer *kstart, integer *kstop, 
	integer *n, real *d__, real *q, integer *ldq, real *rho, real *dlamda,
	 real *w, real *s, integer *lds, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(slaeda)(integer *n, integer *tlvls, integer *curlvl, 
	integer *curpbm, integer *prmptr, integer *perm, integer *givptr, 
	integer *givcol, real *givnum, real *q, integer *qptr, real *z__, 
	real *ztemp, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(slaein)(logical *rightv, logical *noinit, integer *n, 
	real *h__, integer *ldh, real *wr, real *wi, real *vr, real *vi, real 
	*b, integer *ldb, real *work, real *eps3, real *smlnum, real *bignum, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(slaev2)(real *a, real *b, real *c__, real *rt1, real *
	rt2, real *cs1, real *sn1);
 
/* Subroutine */ int LAPACKMANGLE(slaexc)(logical *wantq, integer *n, real *t, integer *
	ldt, real *q, integer *ldq, integer *j1, integer *n1, integer *n2, 
	real *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(slag2)(real *a, integer *lda, real *b, integer *ldb, 
	real *safmin, real *scale1, real *scale2, real *wr1, real *wr2, real *
	wi);
 
/* Subroutine */ int LAPACKMANGLE(slags2)(logical *upper, real *a1, real *a2, real *a3, 
	real *b1, real *b2, real *b3, real *csu, real *snu, real *csv, real *
	snv, real *csq, real *snq);
 
/* Subroutine */ int LAPACKMANGLE(slagtf)(integer *n, real *a, real *lambda, real *b, real 
	*c__, real *tol, real *d__, integer *in, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(slagtm)(char *trans, integer *n, integer *nrhs, real *
	alpha, real *dl, real *d__, real *du, real *x, integer *ldx, real *
	beta, real *b, integer *ldb);
 
/* Subroutine */ int LAPACKMANGLE(slagts)(integer *job, integer *n, real *a, real *b, real 
	*c__, real *d__, integer *in, real *y, real *tol, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(slagv2)(real *a, integer *lda, real *b, integer *ldb, 
	real *alphar, real *alphai, real *beta, real *csl, real *snl, real *
	csr, real *snr);
 
/* Subroutine */ int LAPACKMANGLE(slahqr)(logical *wantt, logical *wantz, integer *n, 
	integer *ilo, integer *ihi, real *h__, integer *ldh, real *wr, real *
	wi, integer *iloz, integer *ihiz, real *z__, integer *ldz, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(slahrd)(integer *n, integer *k, integer *nb, real *a, 
	integer *lda, real *tau, real *t, integer *ldt, real *y, integer *ldy);
 
/* Subroutine */ int LAPACKMANGLE(slaic1)(integer *job, integer *j, real *x, real *sest, 
	real *w, real *gamma, real *sestpr, real *s, real *c__);
 
/* Subroutine */ int LAPACKMANGLE(slaln2)(logical *ltrans, integer *na, integer *nw, real *
	smin, real *ca, real *a, integer *lda, real *d1, real *d2, real *b, 
	integer *ldb, real *wr, real *wi, real *x, integer *ldx, real *scale, 
	real *xnorm, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(slals0)(integer *icompq, integer *nl, integer *nr, 
	integer *sqre, integer *nrhs, real *b, integer *ldb, real *bx, 
	integer *ldbx, integer *perm, integer *givptr, integer *givcol, 
	integer *ldgcol, real *givnum, integer *ldgnum, real *poles, real *
	difl, real *difr, real *z__, integer *k, real *c__, real *s, real *
	work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(slalsa)(integer *icompq, integer *smlsiz, integer *n, 
	integer *nrhs, real *b, integer *ldb, real *bx, integer *ldbx, real *
	u, integer *ldu, real *vt, integer *k, real *difl, real *difr, real *
	z__, real *poles, integer *givptr, integer *givcol, integer *ldgcol, 
	integer *perm, real *givnum, real *c__, real *s, real *work, integer *
	iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(slalsd)(char *uplo, integer *smlsiz, integer *n, integer 
	*nrhs, real *d__, real *e, real *b, integer *ldb, real *rcond, 
	integer *rank, real *work, integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(slamc1)(integer *beta, integer *t, logical *rnd, logical 
	*ieee1);
 
/* Subroutine */ int LAPACKMANGLE(slamc2)(integer *beta, integer *t, logical *rnd, real *
	eps, integer *emin, real *rmin, integer *emax, real *rmax);
 
/* Subroutine */ int LAPACKMANGLE(slamc4)(integer *emin, real *start, integer *base);
 
/* Subroutine */ int LAPACKMANGLE(slamc5)(integer *beta, integer *p, integer *emin, 
	logical *ieee, integer *emax, real *rmax);
 
/* Subroutine */ int LAPACKMANGLE(slamrg)(integer *n1, integer *n2, real *a, integer *
	strd1, integer *strd2, integer *index);
 
/* Subroutine */ int LAPACKMANGLE(slanv2)(real *a, real *b, real *c__, real *d__, real *
	rt1r, real *rt1i, real *rt2r, real *rt2i, real *cs, real *sn);
 
/* Subroutine */ int LAPACKMANGLE(slapll)(integer *n, real *x, integer *incx, real *y, 
	integer *incy, real *ssmin);
 
/* Subroutine */ int LAPACKMANGLE(slapmt)(logical *forwrd, integer *m, integer *n, real *x,
	 integer *ldx, integer *k);
 
/* Subroutine */ int LAPACKMANGLE(slaqgb)(integer *m, integer *n, integer *kl, integer *ku,
	 real *ab, integer *ldab, real *r__, real *c__, real *rowcnd, real *
	colcnd, real *amax, char *equed);
 
/* Subroutine */ int LAPACKMANGLE(slaqge)(integer *m, integer *n, real *a, integer *lda, 
	real *r__, real *c__, real *rowcnd, real *colcnd, real *amax, char *
	equed);
 
/* Subroutine */ int LAPACKMANGLE(slaqp2)(integer *m, integer *n, integer *offset, real *a,
	 integer *lda, integer *jpvt, real *tau, real *vn1, real *vn2, real *
	work);
 
/* Subroutine */ int LAPACKMANGLE(slaqps)(integer *m, integer *n, integer *offset, integer 
	*nb, integer *kb, real *a, integer *lda, integer *jpvt, real *tau, 
	real *vn1, real *vn2, real *auxv, real *f, integer *ldf);
 
/* Subroutine */ int LAPACKMANGLE(slaqsb)(char *uplo, integer *n, integer *kd, real *ab, 
	integer *ldab, real *s, real *scond, real *amax, char *equed);
 
/* Subroutine */ int LAPACKMANGLE(slaqsp)(char *uplo, integer *n, real *ap, real *s, real *
	scond, real *amax, char *equed);
 
/* Subroutine */ int LAPACKMANGLE(slaqsy)(char *uplo, integer *n, real *a, integer *lda, 
	real *s, real *scond, real *amax, char *equed);
 
/* Subroutine */ int LAPACKMANGLE(slaqtr)(logical *ltran, logical *lreal, integer *n, real 
	*t, integer *ldt, real *b, real *w, real *scale, real *x, real *work, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(slar1v)(integer *n, integer *b1, integer *bn, real *
	sigma, real *d__, real *l, real *ld, real *lld, real *gersch, real *
	z__, real *ztz, real *mingma, integer *r__, integer *isuppz, real *
	work);
 
/* Subroutine */ int LAPACKMANGLE(slar2v)(integer *n, real *x, real *y, real *z__, integer 
	*incx, real *c__, real *s, integer *incc);
 
/* Subroutine */ int LAPACKMANGLE(slarf)(char *side, integer *m, integer *n, real *v, 
	integer *incv, real *tau, real *c__, integer *ldc, real *work);
 
/* Subroutine */ int LAPACKMANGLE(slarfb)(char *side, char *trans, char *direct, char *
	storev, integer *m, integer *n, integer *k, real *v, integer *ldv, 
	real *t, integer *ldt, real *c__, integer *ldc, real *work, integer *
	ldwork);
 
/* Subroutine */ int LAPACKMANGLE(slarfg)(integer *n, real *alpha, real *x, integer *incx, 
	real *tau);
 
/* Subroutine */ int LAPACKMANGLE(slarft)(char *direct, char *storev, integer *n, integer *
	k, real *v, integer *ldv, real *tau, real *t, integer *ldt);
 
/* Subroutine */ int LAPACKMANGLE(slarfx)(char *side, integer *m, integer *n, real *v, 
	real *tau, real *c__, integer *ldc, real *work);
 
/* Subroutine */ int LAPACKMANGLE(slargv)(integer *n, real *x, integer *incx, real *y, 
	integer *incy, real *c__, integer *incc);
 
/* Subroutine */ int LAPACKMANGLE(slarnv)(integer *idist, integer *iseed, integer *n, real 
	*x);
 
/* Subroutine */ int LAPACKMANGLE(slarrb)(integer *n, real *d__, real *l, real *ld, real *
	lld, integer *ifirst, integer *ilast, real *sigma, real *reltol, real 
	*w, real *wgap, real *werr, real *work, integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(slarre)(integer *n, real *d__, real *e, real *tol, 
	integer *nsplit, integer *isplit, integer *m, real *w, real *woff, 
	real *gersch, real *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(slarrf)(integer *n, real *d__, real *l, real *ld, real *
	lld, integer *ifirst, integer *ilast, real *w, real *dplus, real *
	lplus, real *work, integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(slarrv)(integer *n, real *d__, real *l, integer *isplit, 
	integer *m, real *w, integer *iblock, real *gersch, real *tol, real *
	z__, integer *ldz, integer *isuppz, real *work, integer *iwork, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(slartg)(real *f, real *g, real *cs, real *sn, real *r__);
 
/* Subroutine */ int LAPACKMANGLE(slartv)(integer *n, real *x, integer *incx, real *y, 
	integer *incy, real *c__, real *s, integer *incc);
 
/* Subroutine */ int LAPACKMANGLE(slaruv)(integer *iseed, integer *n, real *x);
 
/* Subroutine */ int LAPACKMANGLE(slarz)(char *side, integer *m, integer *n, integer *l, 
	real *v, integer *incv, real *tau, real *c__, integer *ldc, real *
	work);
 
/* Subroutine */ int LAPACKMANGLE(slarzb)(char *side, char *trans, char *direct, char *
	storev, integer *m, integer *n, integer *k, integer *l, real *v, 
	integer *ldv, real *t, integer *ldt, real *c__, integer *ldc, real *
	work, integer *ldwork);
 
/* Subroutine */ int LAPACKMANGLE(slarzt)(char *direct, char *storev, integer *n, integer *
	k, real *v, integer *ldv, real *tau, real *t, integer *ldt);
 
/* Subroutine */ int LAPACKMANGLE(slas2)(real *f, real *g, real *h__, real *ssmin, real *
	ssmax);
 
/* Subroutine */ int LAPACKMANGLE(slascl)(char *type__, integer *kl, integer *ku, real *
	cfrom, real *cto, integer *m, integer *n, real *a, integer *lda, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(slasd0)(integer *n, integer *sqre, real *d__, real *e, 
	real *u, integer *ldu, real *vt, integer *ldvt, integer *smlsiz, 
	integer *iwork, real *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(slasd1)(integer *nl, integer *nr, integer *sqre, real *
	d__, real *alpha, real *beta, real *u, integer *ldu, real *vt, 
	integer *ldvt, integer *idxq, integer *iwork, real *work, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(slasd2)(integer *nl, integer *nr, integer *sqre, integer 
	*k, real *d__, real *z__, real *alpha, real *beta, real *u, integer *
	ldu, real *vt, integer *ldvt, real *dsigma, real *u2, integer *ldu2, 
	real *vt2, integer *ldvt2, integer *idxp, integer *idx, integer *idxc,
	 integer *idxq, integer *coltyp, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(slasd3)(integer *nl, integer *nr, integer *sqre, integer 
	*k, real *d__, real *q, integer *ldq, real *dsigma, real *u, integer *
	ldu, real *u2, integer *ldu2, real *vt, integer *ldvt, real *vt2, 
	integer *ldvt2, integer *idxc, integer *ctot, real *z__, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(slasd4)(integer *n, integer *i__, real *d__, real *z__, 
	real *delta, real *rho, real *sigma, real *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(slasd5)(integer *i__, real *d__, real *z__, real *delta, 
	real *rho, real *dsigma, real *work);
 
/* Subroutine */ int LAPACKMANGLE(slasd6)(integer *icompq, integer *nl, integer *nr, 
	integer *sqre, real *d__, real *vf, real *vl, real *alpha, real *beta,
	 integer *idxq, integer *perm, integer *givptr, integer *givcol, 
	integer *ldgcol, real *givnum, integer *ldgnum, real *poles, real *
	difl, real *difr, real *z__, integer *k, real *c__, real *s, real *
	work, integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(slasd7)(integer *icompq, integer *nl, integer *nr, 
	integer *sqre, integer *k, real *d__, real *z__, real *zw, real *vf, 
	real *vfw, real *vl, real *vlw, real *alpha, real *beta, real *dsigma,
	 integer *idx, integer *idxp, integer *idxq, integer *perm, integer *
	givptr, integer *givcol, integer *ldgcol, real *givnum, integer *
	ldgnum, real *c__, real *s, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(slasd8)(integer *icompq, integer *k, real *d__, real *
	z__, real *vf, real *vl, real *difl, real *difr, integer *lddifr, 
	real *dsigma, real *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(slasd9)(integer *icompq, integer *ldu, integer *k, real *
	d__, real *z__, real *vf, real *vl, real *difl, real *difr, real *
	dsigma, real *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(slasda)(integer *icompq, integer *smlsiz, integer *n, 
	integer *sqre, real *d__, real *e, real *u, integer *ldu, real *vt, 
	integer *k, real *difl, real *difr, real *z__, real *poles, integer *
	givptr, integer *givcol, integer *ldgcol, integer *perm, real *givnum,
	 real *c__, real *s, real *work, integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(slasdq)(char *uplo, integer *sqre, integer *n, integer *
	ncvt, integer *nru, integer *ncc, real *d__, real *e, real *vt, 
	integer *ldvt, real *u, integer *ldu, real *c__, integer *ldc, real *
	work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(slasdt)(integer *n, integer *lvl, integer *nd, integer *
	inode, integer *ndiml, integer *ndimr, integer *msub);
 
/* Subroutine */ int LAPACKMANGLE(slaset)(char *uplo, integer *m, integer *n, real *alpha, 
	real *beta, real *a, integer *lda);
 
/* Subroutine */ int LAPACKMANGLE(slasq1)(integer *n, real *d__, real *e, real *work, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(slasq2)(integer *n, real *z__, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(slasq3)(integer *i0, integer *n0, real *z__, integer *pp,
	 real *dmin__, real *sigma, real *desig, real *qmax, integer *nfail, 
	integer *iter, integer *ndiv, logical *ieee);
 
/* Subroutine */ int LAPACKMANGLE(slasq4)(integer *i0, integer *n0, real *z__, integer *pp,
	 integer *n0in, real *dmin__, real *dmin1, real *dmin2, real *dn, 
	real *dn1, real *dn2, real *tau, integer *ttype);
 
/* Subroutine */ int LAPACKMANGLE(slasq5)(integer *i0, integer *n0, real *z__, integer *pp,
	 real *tau, real *dmin__, real *dmin1, real *dmin2, real *dn, real *
	dnm1, real *dnm2, logical *ieee);
 
/* Subroutine */ int LAPACKMANGLE(slasq6)(integer *i0, integer *n0, real *z__, integer *pp,
	 real *dmin__, real *dmin1, real *dmin2, real *dn, real *dnm1, real *
	dnm2);
 
/* Subroutine */ int LAPACKMANGLE(slasr)(char *side, char *pivot, char *direct, integer *m,
	 integer *n, real *c__, real *s, real *a, integer *lda);
 
/* Subroutine */ int LAPACKMANGLE(slasrt)(char *id, integer *n, real *d__, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(slassq)(integer *n, real *x, integer *incx, real *scale, 
	real *sumsq);
 
/* Subroutine */ int LAPACKMANGLE(slasv2)(real *f, real *g, real *h__, real *ssmin, real *
	ssmax, real *snr, real *csr, real *snl, real *csl);
 
/* Subroutine */ int LAPACKMANGLE(slaswp)(integer *n, real *a, integer *lda, integer *k1, 
	integer *k2, integer *ipiv, integer *incx);
 
/* Subroutine */ int LAPACKMANGLE(slasy2)(logical *ltranl, logical *ltranr, integer *isgn, 
	integer *n1, integer *n2, real *tl, integer *ldtl, real *tr, integer *
	ldtr, real *b, integer *ldb, real *scale, real *x, integer *ldx, real 
	*xnorm, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(slasyf)(char *uplo, integer *n, integer *nb, integer *kb,
	 real *a, integer *lda, integer *ipiv, real *w, integer *ldw, integer 
	*info);
 
/* Subroutine */ int LAPACKMANGLE(slatbs)(char *uplo, char *trans, char *diag, char *
	normin, integer *n, integer *kd, real *ab, integer *ldab, real *x, 
	real *scale, real *cnorm, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(slatdf)(integer *ijob, integer *n, real *z__, integer *
	ldz, real *rhs, real *rdsum, real *rdscal, integer *ipiv, integer *
	jpiv);
 
/* Subroutine */ int LAPACKMANGLE(slatps)(char *uplo, char *trans, char *diag, char *
	normin, integer *n, real *ap, real *x, real *scale, real *cnorm, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(slatrd)(char *uplo, integer *n, integer *nb, real *a, 
	integer *lda, real *e, real *tau, real *w, integer *ldw);
 
/* Subroutine */ int LAPACKMANGLE(slatrs)(char *uplo, char *trans, char *diag, char *
	normin, integer *n, real *a, integer *lda, real *x, real *scale, real 
	*cnorm, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(slatrz)(integer *m, integer *n, integer *l, real *a, 
	integer *lda, real *tau, real *work);
 
/* Subroutine */ int LAPACKMANGLE(slatzm)(char *side, integer *m, integer *n, real *v, 
	integer *incv, real *tau, real *c1, real *c2, integer *ldc, real *
	work);
 
/* Subroutine */ int LAPACKMANGLE(slauu2)(char *uplo, integer *n, real *a, integer *lda, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(slauum)(char *uplo, integer *n, real *a, integer *lda, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sopgtr)(char *uplo, integer *n, real *ap, real *tau, 
	real *q, integer *ldq, real *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sopmtr)(char *side, char *uplo, char *trans, integer *m, 
	integer *n, real *ap, real *tau, real *c__, integer *ldc, real *work, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sorg2l)(integer *m, integer *n, integer *k, real *a, 
	integer *lda, real *tau, real *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sorg2r)(integer *m, integer *n, integer *k, real *a, 
	integer *lda, real *tau, real *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sorgbr)(char *vect, integer *m, integer *n, integer *k, 
	real *a, integer *lda, real *tau, real *work, integer *lwork, integer 
	*info);
 
/* Subroutine */ int LAPACKMANGLE(sorghr)(integer *n, integer *ilo, integer *ihi, real *a, 
	integer *lda, real *tau, real *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sorgl2)(integer *m, integer *n, integer *k, real *a, 
	integer *lda, real *tau, real *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sorglq)(integer *m, integer *n, integer *k, real *a, 
	integer *lda, real *tau, real *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sorgql)(integer *m, integer *n, integer *k, real *a, 
	integer *lda, real *tau, real *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sorgqr)(integer *m, integer *n, integer *k, real *a, 
	integer *lda, real *tau, real *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sorgr2)(integer *m, integer *n, integer *k, real *a, 
	integer *lda, real *tau, real *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sorgrq)(integer *m, integer *n, integer *k, real *a, 
	integer *lda, real *tau, real *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sorgtr)(char *uplo, integer *n, real *a, integer *lda, 
	real *tau, real *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sorm2l)(char *side, char *trans, integer *m, integer *n, 
	integer *k, real *a, integer *lda, real *tau, real *c__, integer *ldc,
	 real *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sorm2r)(char *side, char *trans, integer *m, integer *n, 
	integer *k, real *a, integer *lda, real *tau, real *c__, integer *ldc,
	 real *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sormbr)(char *vect, char *side, char *trans, integer *m, 
	integer *n, integer *k, real *a, integer *lda, real *tau, real *c__, 
	integer *ldc, real *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sormhr)(char *side, char *trans, integer *m, integer *n, 
	integer *ilo, integer *ihi, real *a, integer *lda, real *tau, real *
	c__, integer *ldc, real *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sorml2)(char *side, char *trans, integer *m, integer *n, 
	integer *k, real *a, integer *lda, real *tau, real *c__, integer *ldc,
	 real *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sormlq)(char *side, char *trans, integer *m, integer *n, 
	integer *k, real *a, integer *lda, real *tau, real *c__, integer *ldc,
	 real *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sormql)(char *side, char *trans, integer *m, integer *n, 
	integer *k, real *a, integer *lda, real *tau, real *c__, integer *ldc,
	 real *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sormqr)(char *side, char *trans, integer *m, integer *n, 
	integer *k, real *a, integer *lda, real *tau, real *c__, integer *ldc,
	 real *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sormr2)(char *side, char *trans, integer *m, integer *n, 
	integer *k, real *a, integer *lda, real *tau, real *c__, integer *ldc,
	 real *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sormr3)(char *side, char *trans, integer *m, integer *n, 
	integer *k, integer *l, real *a, integer *lda, real *tau, real *c__, 
	integer *ldc, real *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sormrq)(char *side, char *trans, integer *m, integer *n, 
	integer *k, real *a, integer *lda, real *tau, real *c__, integer *ldc,
	 real *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sormrz)(char *side, char *trans, integer *m, integer *n, 
	integer *k, integer *l, real *a, integer *lda, real *tau, real *c__, 
	integer *ldc, real *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sormtr)(char *side, char *uplo, char *trans, integer *m, 
	integer *n, real *a, integer *lda, real *tau, real *c__, integer *ldc,
	 real *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(spbcon)(char *uplo, integer *n, integer *kd, real *ab, 
	integer *ldab, real *anorm, real *rcond, real *work, integer *iwork, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(spbequ)(char *uplo, integer *n, integer *kd, real *ab, 
	integer *ldab, real *s, real *scond, real *amax, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(spbrfs)(char *uplo, integer *n, integer *kd, integer *
	nrhs, real *ab, integer *ldab, real *afb, integer *ldafb, real *b, 
	integer *ldb, real *x, integer *ldx, real *ferr, real *berr, real *
	work, integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(spbstf)(char *uplo, integer *n, integer *kd, real *ab, 
	integer *ldab, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(spbsv)(char *uplo, integer *n, integer *kd, integer *
	nrhs, real *ab, integer *ldab, real *b, integer *ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(spbsvx)(char *fact, char *uplo, integer *n, integer *kd, 
	integer *nrhs, real *ab, integer *ldab, real *afb, integer *ldafb, 
	char *equed, real *s, real *b, integer *ldb, real *x, integer *ldx, 
	real *rcond, real *ferr, real *berr, real *work, integer *iwork, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(spbtf2)(char *uplo, integer *n, integer *kd, real *ab, 
	integer *ldab, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(spbtrf)(char *uplo, integer *n, integer *kd, real *ab, 
	integer *ldab, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(spbtrs)(char *uplo, integer *n, integer *kd, integer *
	nrhs, real *ab, integer *ldab, real *b, integer *ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(spocon)(char *uplo, integer *n, real *a, integer *lda, 
	real *anorm, real *rcond, real *work, integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(spoequ)(integer *n, real *a, integer *lda, real *s, real 
	*scond, real *amax, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sporfs)(char *uplo, integer *n, integer *nrhs, real *a, 
	integer *lda, real *af, integer *ldaf, real *b, integer *ldb, real *x,
	 integer *ldx, real *ferr, real *berr, real *work, integer *iwork, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sposv)(char *uplo, integer *n, integer *nrhs, real *a, 
	integer *lda, real *b, integer *ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sposvx)(char *fact, char *uplo, integer *n, integer *
	nrhs, real *a, integer *lda, real *af, integer *ldaf, char *equed, 
	real *s, real *b, integer *ldb, real *x, integer *ldx, real *rcond, 
	real *ferr, real *berr, real *work, integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(spotf2)(char *uplo, integer *n, real *a, integer *lda, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(spotrf)(char *uplo, integer *n, real *a, integer *lda, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(spotri)(char *uplo, integer *n, real *a, integer *lda, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(spotrs)(char *uplo, integer *n, integer *nrhs, real *a, 
	integer *lda, real *b, integer *ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sppcon)(char *uplo, integer *n, real *ap, real *anorm, 
	real *rcond, real *work, integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sppequ)(char *uplo, integer *n, real *ap, real *s, real *
	scond, real *amax, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(spprfs)(char *uplo, integer *n, integer *nrhs, real *ap, 
	real *afp, real *b, integer *ldb, real *x, integer *ldx, real *ferr, 
	real *berr, real *work, integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sppsv)(char *uplo, integer *n, integer *nrhs, real *ap, 
	real *b, integer *ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sppsvx)(char *fact, char *uplo, integer *n, integer *
	nrhs, real *ap, real *afp, char *equed, real *s, real *b, integer *
	ldb, real *x, integer *ldx, real *rcond, real *ferr, real *berr, real 
	*work, integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(spptrf)(char *uplo, integer *n, real *ap, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(spptri)(char *uplo, integer *n, real *ap, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(spptrs)(char *uplo, integer *n, integer *nrhs, real *ap, 
	real *b, integer *ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sptcon)(integer *n, real *d__, real *e, real *anorm, 
	real *rcond, real *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(spteqr)(char *compz, integer *n, real *d__, real *e, 
	real *z__, integer *ldz, real *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sptrfs)(integer *n, integer *nrhs, real *d__, real *e, 
	real *df, real *ef, real *b, integer *ldb, real *x, integer *ldx, 
	real *ferr, real *berr, real *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sptsv)(integer *n, integer *nrhs, real *d__, real *e, 
	real *b, integer *ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sptsvx)(char *fact, integer *n, integer *nrhs, real *d__,
	 real *e, real *df, real *ef, real *b, integer *ldb, real *x, integer 
	*ldx, real *rcond, real *ferr, real *berr, real *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(spttrf)(integer *n, real *d__, real *e, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(spttrs)(integer *n, integer *nrhs, real *d__, real *e, 
	real *b, integer *ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sptts2)(integer *n, integer *nrhs, real *d__, real *e, 
	real *b, integer *ldb);
 
/* Subroutine */ int LAPACKMANGLE(srscl)(integer *n, real *sa, real *sx, integer *incx);
 
/* Subroutine */ int LAPACKMANGLE(ssbev)(char *jobz, char *uplo, integer *n, integer *kd, 
	real *ab, integer *ldab, real *w, real *z__, integer *ldz, real *work,
	 integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ssbevd)(char *jobz, char *uplo, integer *n, integer *kd, 
	real *ab, integer *ldab, real *w, real *z__, integer *ldz, real *work,
	 integer *lwork, integer *iwork, integer *liwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ssbevx)(char *jobz, char *range, char *uplo, integer *n, 
	integer *kd, real *ab, integer *ldab, real *q, integer *ldq, real *vl,
	 real *vu, integer *il, integer *iu, real *abstol, integer *m, real *
	w, real *z__, integer *ldz, real *work, integer *iwork, integer *
	ifail, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ssbgst)(char *vect, char *uplo, integer *n, integer *ka, 
	integer *kb, real *ab, integer *ldab, real *bb, integer *ldbb, real *
	x, integer *ldx, real *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ssbgv)(char *jobz, char *uplo, integer *n, integer *ka, 
	integer *kb, real *ab, integer *ldab, real *bb, integer *ldbb, real *
	w, real *z__, integer *ldz, real *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ssbgvd)(char *jobz, char *uplo, integer *n, integer *ka, 
	integer *kb, real *ab, integer *ldab, real *bb, integer *ldbb, real *
	w, real *z__, integer *ldz, real *work, integer *lwork, integer *
	iwork, integer *liwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ssbgvx)(char *jobz, char *range, char *uplo, integer *n, 
	integer *ka, integer *kb, real *ab, integer *ldab, real *bb, integer *
	ldbb, real *q, integer *ldq, real *vl, real *vu, integer *il, integer 
	*iu, real *abstol, integer *m, real *w, real *z__, integer *ldz, real 
	*work, integer *iwork, integer *ifail, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ssbtrd)(char *vect, char *uplo, integer *n, integer *kd, 
	real *ab, integer *ldab, real *d__, real *e, real *q, integer *ldq, 
	real *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sspcon)(char *uplo, integer *n, real *ap, integer *ipiv, 
	real *anorm, real *rcond, real *work, integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sspev)(char *jobz, char *uplo, integer *n, real *ap, 
	real *w, real *z__, integer *ldz, real *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sspevd)(char *jobz, char *uplo, integer *n, real *ap, 
	real *w, real *z__, integer *ldz, real *work, integer *lwork, integer 
	*iwork, integer *liwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sspevx)(char *jobz, char *range, char *uplo, integer *n, 
	real *ap, real *vl, real *vu, integer *il, integer *iu, real *abstol, 
	integer *m, real *w, real *z__, integer *ldz, real *work, integer *
	iwork, integer *ifail, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sspgst)(integer *itype, char *uplo, integer *n, real *ap,
	 real *bp, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sspgv)(integer *itype, char *jobz, char *uplo, integer *
	n, real *ap, real *bp, real *w, real *z__, integer *ldz, real *work, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sspgvd)(integer *itype, char *jobz, char *uplo, integer *
	n, real *ap, real *bp, real *w, real *z__, integer *ldz, real *work, 
	integer *lwork, integer *iwork, integer *liwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sspgvx)(integer *itype, char *jobz, char *range, char *
	uplo, integer *n, real *ap, real *bp, real *vl, real *vu, integer *il,
	 integer *iu, real *abstol, integer *m, real *w, real *z__, integer *
	ldz, real *work, integer *iwork, integer *ifail, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ssprfs)(char *uplo, integer *n, integer *nrhs, real *ap, 
	real *afp, integer *ipiv, real *b, integer *ldb, real *x, integer *
	ldx, real *ferr, real *berr, real *work, integer *iwork, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(sspsv)(char *uplo, integer *n, integer *nrhs, real *ap, 
	integer *ipiv, real *b, integer *ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sspsvx)(char *fact, char *uplo, integer *n, integer *
	nrhs, real *ap, real *afp, integer *ipiv, real *b, integer *ldb, real 
	*x, integer *ldx, real *rcond, real *ferr, real *berr, real *work, 
	integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ssptrd)(char *uplo, integer *n, real *ap, real *d__, 
	real *e, real *tau, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ssptrf)(char *uplo, integer *n, real *ap, integer *ipiv, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ssptri)(char *uplo, integer *n, real *ap, integer *ipiv, 
	real *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ssptrs)(char *uplo, integer *n, integer *nrhs, real *ap, 
	integer *ipiv, real *b, integer *ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sstebz)(char *range, char *order, integer *n, real *vl, 
	real *vu, integer *il, integer *iu, real *abstol, real *d__, real *e, 
	integer *m, integer *nsplit, real *w, integer *iblock, integer *
	isplit, real *work, integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sstedc)(char *compz, integer *n, real *d__, real *e, 
	real *z__, integer *ldz, real *work, integer *lwork, integer *iwork, 
	integer *liwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sstegr)(char *jobz, char *range, integer *n, real *d__, 
	real *e, real *vl, real *vu, integer *il, integer *iu, real *abstol, 
	integer *m, real *w, real *z__, integer *ldz, integer *isuppz, real *
	work, integer *lwork, integer *iwork, integer *liwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sstein)(integer *n, real *d__, real *e, integer *m, real 
	*w, integer *iblock, integer *isplit, real *z__, integer *ldz, real *
	work, integer *iwork, integer *ifail, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ssteqr)(char *compz, integer *n, real *d__, real *e, 
	real *z__, integer *ldz, real *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ssterf)(integer *n, real *d__, real *e, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sstev)(char *jobz, integer *n, real *d__, real *e, real *
	z__, integer *ldz, real *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sstevd)(char *jobz, integer *n, real *d__, real *e, real 
	*z__, integer *ldz, real *work, integer *lwork, integer *iwork, 
	integer *liwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sstevr)(char *jobz, char *range, integer *n, real *d__, 
	real *e, real *vl, real *vu, integer *il, integer *iu, real *abstol, 
	integer *m, real *w, real *z__, integer *ldz, integer *isuppz, real *
	work, integer *lwork, integer *iwork, integer *liwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(sstevx)(char *jobz, char *range, integer *n, real *d__, 
	real *e, real *vl, real *vu, integer *il, integer *iu, real *abstol, 
	integer *m, real *w, real *z__, integer *ldz, real *work, integer *
	iwork, integer *ifail, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ssycon)(char *uplo, integer *n, real *a, integer *lda, 
	integer *ipiv, real *anorm, real *rcond, real *work, integer *iwork, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ssyev)(char *jobz, char *uplo, integer *n, real *a, 
	integer *lda, real *w, real *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ssyevd)(char *jobz, char *uplo, integer *n, real *a, 
	integer *lda, real *w, real *work, integer *lwork, integer *iwork, 
	integer *liwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ssyevr)(char *jobz, char *range, char *uplo, integer *n, 
	real *a, integer *lda, real *vl, real *vu, integer *il, integer *iu, 
	real *abstol, integer *m, real *w, real *z__, integer *ldz, integer *
	isuppz, real *work, integer *lwork, integer *iwork, integer *liwork, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ssyevx)(char *jobz, char *range, char *uplo, integer *n, 
	real *a, integer *lda, real *vl, real *vu, integer *il, integer *iu, 
	real *abstol, integer *m, real *w, real *z__, integer *ldz, real *
	work, integer *lwork, integer *iwork, integer *ifail, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ssygs2)(integer *itype, char *uplo, integer *n, real *a, 
	integer *lda, real *b, integer *ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ssygst)(integer *itype, char *uplo, integer *n, real *a, 
	integer *lda, real *b, integer *ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ssygv)(integer *itype, char *jobz, char *uplo, integer *
	n, real *a, integer *lda, real *b, integer *ldb, real *w, real *work, 
	integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ssygvd)(integer *itype, char *jobz, char *uplo, integer *
	n, real *a, integer *lda, real *b, integer *ldb, real *w, real *work, 
	integer *lwork, integer *iwork, integer *liwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ssygvx)(integer *itype, char *jobz, char *range, char *
	uplo, integer *n, real *a, integer *lda, real *b, integer *ldb, real *
	vl, real *vu, integer *il, integer *iu, real *abstol, integer *m, 
	real *w, real *z__, integer *ldz, real *work, integer *lwork, integer 
	*iwork, integer *ifail, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ssyrfs)(char *uplo, integer *n, integer *nrhs, real *a, 
	integer *lda, real *af, integer *ldaf, integer *ipiv, real *b, 
	integer *ldb, real *x, integer *ldx, real *ferr, real *berr, real *
	work, integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ssysv)(char *uplo, integer *n, integer *nrhs, real *a, 
	integer *lda, integer *ipiv, real *b, integer *ldb, real *work, 
	integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ssysvx)(char *fact, char *uplo, integer *n, integer *
	nrhs, real *a, integer *lda, real *af, integer *ldaf, integer *ipiv, 
	real *b, integer *ldb, real *x, integer *ldx, real *rcond, real *ferr,
	 real *berr, real *work, integer *lwork, integer *iwork, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(ssytd2)(char *uplo, integer *n, real *a, integer *lda, 
	real *d__, real *e, real *tau, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ssytf2)(char *uplo, integer *n, real *a, integer *lda, 
	integer *ipiv, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ssytrd)(char *uplo, integer *n, real *a, integer *lda, 
	real *d__, real *e, real *tau, real *work, integer *lwork, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(ssytrf)(char *uplo, integer *n, real *a, integer *lda, 
	integer *ipiv, real *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ssytri)(char *uplo, integer *n, real *a, integer *lda, 
	integer *ipiv, real *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ssytrs)(char *uplo, integer *n, integer *nrhs, real *a, 
	integer *lda, integer *ipiv, real *b, integer *ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(stbcon)(char *norm, char *uplo, char *diag, integer *n, 
	integer *kd, real *ab, integer *ldab, real *rcond, real *work, 
	integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(stbrfs)(char *uplo, char *trans, char *diag, integer *n, 
	integer *kd, integer *nrhs, real *ab, integer *ldab, real *b, integer 
	*ldb, real *x, integer *ldx, real *ferr, real *berr, real *work, 
	integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(stbtrs)(char *uplo, char *trans, char *diag, integer *n, 
	integer *kd, integer *nrhs, real *ab, integer *ldab, real *b, integer 
	*ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(stgevc)(char *side, char *howmny, logical *select, 
	integer *n, real *a, integer *lda, real *b, integer *ldb, real *vl, 
	integer *ldvl, real *vr, integer *ldvr, integer *mm, integer *m, real 
	*work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(stgex2)(logical *wantq, logical *wantz, integer *n, real 
	*a, integer *lda, real *b, integer *ldb, real *q, integer *ldq, real *
	z__, integer *ldz, integer *j1, integer *n1, integer *n2, real *work, 
	integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(stgexc)(logical *wantq, logical *wantz, integer *n, real 
	*a, integer *lda, real *b, integer *ldb, real *q, integer *ldq, real *
	z__, integer *ldz, integer *ifst, integer *ilst, real *work, integer *
	lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(stgsen)(integer *ijob, logical *wantq, logical *wantz, 
	logical *select, integer *n, real *a, integer *lda, real *b, integer *
	ldb, real *alphar, real *alphai, real *beta, real *q, integer *ldq, 
	real *z__, integer *ldz, integer *m, real *pl, real *pr, real *dif, 
	real *work, integer *lwork, integer *iwork, integer *liwork, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(stgsja)(char *jobu, char *jobv, char *jobq, integer *m, 
	integer *p, integer *n, integer *k, integer *l, real *a, integer *lda,
	 real *b, integer *ldb, real *tola, real *tolb, real *alpha, real *
	beta, real *u, integer *ldu, real *v, integer *ldv, real *q, integer *
	ldq, real *work, integer *ncycle, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(stgsna)(char *job, char *howmny, logical *select, 
	integer *n, real *a, integer *lda, real *b, integer *ldb, real *vl, 
	integer *ldvl, real *vr, integer *ldvr, real *s, real *dif, integer *
	mm, integer *m, real *work, integer *lwork, integer *iwork, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(stgsy2)(char *trans, integer *ijob, integer *m, integer *
	n, real *a, integer *lda, real *b, integer *ldb, real *c__, integer *
	ldc, real *d__, integer *ldd, real *e, integer *lde, real *f, integer 
	*ldf, real *scale, real *rdsum, real *rdscal, integer *iwork, integer 
	*pq, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(stgsyl)(char *trans, integer *ijob, integer *m, integer *
	n, real *a, integer *lda, real *b, integer *ldb, real *c__, integer *
	ldc, real *d__, integer *ldd, real *e, integer *lde, real *f, integer 
	*ldf, real *scale, real *dif, real *work, integer *lwork, integer *
	iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(stpcon)(char *norm, char *uplo, char *diag, integer *n, 
	real *ap, real *rcond, real *work, integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(stprfs)(char *uplo, char *trans, char *diag, integer *n, 
	integer *nrhs, real *ap, real *b, integer *ldb, real *x, integer *ldx,
	 real *ferr, real *berr, real *work, integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(stptri)(char *uplo, char *diag, integer *n, real *ap, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(stptrs)(char *uplo, char *trans, char *diag, integer *n, 
	integer *nrhs, real *ap, real *b, integer *ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(strcon)(char *norm, char *uplo, char *diag, integer *n, 
	real *a, integer *lda, real *rcond, real *work, integer *iwork, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(strevc)(char *side, char *howmny, logical *select, 
	integer *n, real *t, integer *ldt, real *vl, integer *ldvl, real *vr, 
	integer *ldvr, integer *mm, integer *m, real *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(strexc)(char *compq, integer *n, real *t, integer *ldt, 
	real *q, integer *ldq, integer *ifst, integer *ilst, real *work, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(strrfs)(char *uplo, char *trans, char *diag, integer *n, 
	integer *nrhs, real *a, integer *lda, real *b, integer *ldb, real *x, 
	integer *ldx, real *ferr, real *berr, real *work, integer *iwork, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(strsen)(char *job, char *compq, logical *select, integer 
	*n, real *t, integer *ldt, real *q, integer *ldq, real *wr, real *wi, 
	integer *m, real *s, real *sep, real *work, integer *lwork, integer *
	iwork, integer *liwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(strsna)(char *job, char *howmny, logical *select, 
	integer *n, real *t, integer *ldt, real *vl, integer *ldvl, real *vr, 
	integer *ldvr, real *s, real *sep, integer *mm, integer *m, real *
	work, integer *ldwork, integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(strsyl)(char *trana, char *tranb, integer *isgn, integer 
	*m, integer *n, real *a, integer *lda, real *b, integer *ldb, real *
	c__, integer *ldc, real *scale, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(strti2)(char *uplo, char *diag, integer *n, real *a, 
	integer *lda, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(strtri)(char *uplo, char *diag, integer *n, real *a, 
	integer *lda, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(strtrs)(char *uplo, char *trans, char *diag, integer *n, 
	integer *nrhs, real *a, integer *lda, real *b, integer *ldb, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(stzrqf)(integer *m, integer *n, real *a, integer *lda, 
	real *tau, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(stzrzf)(integer *m, integer *n, real *a, integer *lda, 
	real *tau, real *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(xerbla)(char *srname, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zbdsqr)(char *uplo, integer *n, integer *ncvt, integer *
	nru, integer *ncc, doublereal *d__, doublereal *e, doublecomplex *vt, 
	integer *ldvt, doublecomplex *u, integer *ldu, doublecomplex *c__, 
	integer *ldc, doublereal *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zdrot)(integer *n, doublecomplex *cx, integer *incx, 
	doublecomplex *cy, integer *incy, doublereal *c__, doublereal *s);
 
/* Subroutine */ int LAPACKMANGLE(zdrscl)(integer *n, doublereal *sa, doublecomplex *sx, 
	integer *incx);
 
/* Subroutine */ int LAPACKMANGLE(zgbbrd)(char *vect, integer *m, integer *n, integer *ncc,
	 integer *kl, integer *ku, doublecomplex *ab, integer *ldab, 
	doublereal *d__, doublereal *e, doublecomplex *q, integer *ldq, 
	doublecomplex *pt, integer *ldpt, doublecomplex *c__, integer *ldc, 
	doublecomplex *work, doublereal *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zgbcon)(char *norm, integer *n, integer *kl, integer *ku,
	 doublecomplex *ab, integer *ldab, integer *ipiv, doublereal *anorm, 
	doublereal *rcond, doublecomplex *work, doublereal *rwork, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(zgbequ)(integer *m, integer *n, integer *kl, integer *ku,
	 doublecomplex *ab, integer *ldab, doublereal *r__, doublereal *c__, 
	doublereal *rowcnd, doublereal *colcnd, doublereal *amax, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(zgbrfs)(char *trans, integer *n, integer *kl, integer *
	ku, integer *nrhs, doublecomplex *ab, integer *ldab, doublecomplex *
	afb, integer *ldafb, integer *ipiv, doublecomplex *b, integer *ldb, 
	doublecomplex *x, integer *ldx, doublereal *ferr, doublereal *berr, 
	doublecomplex *work, doublereal *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zgbsv)(integer *n, integer *kl, integer *ku, integer *
	nrhs, doublecomplex *ab, integer *ldab, integer *ipiv, doublecomplex *
	b, integer *ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zgbsvx)(char *fact, char *trans, integer *n, integer *kl,
	 integer *ku, integer *nrhs, doublecomplex *ab, integer *ldab, 
	doublecomplex *afb, integer *ldafb, integer *ipiv, char *equed, 
	doublereal *r__, doublereal *c__, doublecomplex *b, integer *ldb, 
	doublecomplex *x, integer *ldx, doublereal *rcond, doublereal *ferr, 
	doublereal *berr, doublecomplex *work, doublereal *rwork, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(zgbtf2)(integer *m, integer *n, integer *kl, integer *ku,
	 doublecomplex *ab, integer *ldab, integer *ipiv, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zgbtrf)(integer *m, integer *n, integer *kl, integer *ku,
	 doublecomplex *ab, integer *ldab, integer *ipiv, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zgbtrs)(char *trans, integer *n, integer *kl, integer *
	ku, integer *nrhs, doublecomplex *ab, integer *ldab, integer *ipiv, 
	doublecomplex *b, integer *ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zgebak)(char *job, char *side, integer *n, integer *ilo, 
	integer *ihi, doublereal *scale, integer *m, doublecomplex *v, 
	integer *ldv, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zgebal)(char *job, integer *n, doublecomplex *a, integer 
	*lda, integer *ilo, integer *ihi, doublereal *scale, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zgebd2)(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublereal *d__, doublereal *e, doublecomplex *tauq, 
	doublecomplex *taup, doublecomplex *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zgebrd)(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublereal *d__, doublereal *e, doublecomplex *tauq, 
	doublecomplex *taup, doublecomplex *work, integer *lwork, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(zgecon)(char *norm, integer *n, doublecomplex *a, 
	integer *lda, doublereal *anorm, doublereal *rcond, doublecomplex *
	work, doublereal *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zgeequ)(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublereal *r__, doublereal *c__, doublereal *rowcnd, 
	doublereal *colcnd, doublereal *amax, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zgees)(char *jobvs, char *sort, L_fp select, integer *n, 
	doublecomplex *a, integer *lda, integer *sdim, doublecomplex *w, 
	doublecomplex *vs, integer *ldvs, doublecomplex *work, integer *lwork,
	 doublereal *rwork, logical *bwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zgeesx)(char *jobvs, char *sort, L_fp select, char *
	sense, integer *n, doublecomplex *a, integer *lda, integer *sdim, 
	doublecomplex *w, doublecomplex *vs, integer *ldvs, doublereal *
	rconde, doublereal *rcondv, doublecomplex *work, integer *lwork, 
	doublereal *rwork, logical *bwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zgeev)(char *jobvl, char *jobvr, integer *n, 
	doublecomplex *a, integer *lda, doublecomplex *w, doublecomplex *vl, 
	integer *ldvl, doublecomplex *vr, integer *ldvr, doublecomplex *work, 
	integer *lwork, doublereal *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zgeevx)(char *balanc, char *jobvl, char *jobvr, char *
	sense, integer *n, doublecomplex *a, integer *lda, doublecomplex *w, 
	doublecomplex *vl, integer *ldvl, doublecomplex *vr, integer *ldvr, 
	integer *ilo, integer *ihi, doublereal *scale, doublereal *abnrm, 
	doublereal *rconde, doublereal *rcondv, doublecomplex *work, integer *
	lwork, doublereal *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zgegs)(char *jobvsl, char *jobvsr, integer *n, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublecomplex *alpha, doublecomplex *beta, doublecomplex *vsl, 
	integer *ldvsl, doublecomplex *vsr, integer *ldvsr, doublecomplex *
	work, integer *lwork, doublereal *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zgegv)(char *jobvl, char *jobvr, integer *n, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublecomplex *alpha, doublecomplex *beta, doublecomplex *vl, integer 
	*ldvl, doublecomplex *vr, integer *ldvr, doublecomplex *work, integer 
	*lwork, doublereal *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zgehd2)(integer *n, integer *ilo, integer *ihi, 
	doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
	work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zgehrd)(integer *n, integer *ilo, integer *ihi, 
	doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
	work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zgelq2)(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublecomplex *tau, doublecomplex *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zgelqf)(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublecomplex *tau, doublecomplex *work, integer *lwork,
	 integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zgels)(char *trans, integer *m, integer *n, integer *
	nrhs, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublecomplex *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zgelsx)(integer *m, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	integer *jpvt, doublereal *rcond, integer *rank, doublecomplex *work, 
	doublereal *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zgelsy)(integer *m, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	integer *jpvt, doublereal *rcond, integer *rank, doublecomplex *work, 
	integer *lwork, doublereal *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zgeql2)(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublecomplex *tau, doublecomplex *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zgeqlf)(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublecomplex *tau, doublecomplex *work, integer *lwork,
	 integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zgeqp3)(integer *m, integer *n, doublecomplex *a, 
	integer *lda, integer *jpvt, doublecomplex *tau, doublecomplex *work, 
	integer *lwork, doublereal *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zgeqpf)(integer *m, integer *n, doublecomplex *a, 
	integer *lda, integer *jpvt, doublecomplex *tau, doublecomplex *work, 
	doublereal *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zgeqr2)(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublecomplex *tau, doublecomplex *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zgeqrf)(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublecomplex *tau, doublecomplex *work, integer *lwork,
	 integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zgerfs)(char *trans, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, doublecomplex *af, integer *ldaf, 
	integer *ipiv, doublecomplex *b, integer *ldb, doublecomplex *x, 
	integer *ldx, doublereal *ferr, doublereal *berr, doublecomplex *work,
	 doublereal *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zgerq2)(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublecomplex *tau, doublecomplex *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zgerqf)(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublecomplex *tau, doublecomplex *work, integer *lwork,
	 integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zgesc2)(integer *n, doublecomplex *a, integer *lda, 
	doublecomplex *rhs, integer *ipiv, integer *jpiv, doublereal *scale);
 
/* Subroutine */ int LAPACKMANGLE(zgesv)(integer *n, integer *nrhs, doublecomplex *a, 
	integer *lda, integer *ipiv, doublecomplex *b, integer *ldb, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(zgesvx)(char *fact, char *trans, integer *n, integer *
	nrhs, doublecomplex *a, integer *lda, doublecomplex *af, integer *
	ldaf, integer *ipiv, char *equed, doublereal *r__, doublereal *c__, 
	doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, 
	doublereal *rcond, doublereal *ferr, doublereal *berr, doublecomplex *
	work, doublereal *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zgetc2)(integer *n, doublecomplex *a, integer *lda, 
	integer *ipiv, integer *jpiv, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zgetf2)(integer *m, integer *n, doublecomplex *a, 
	integer *lda, integer *ipiv, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zgetrf)(integer *m, integer *n, doublecomplex *a, 
	integer *lda, integer *ipiv, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zgetri)(integer *n, doublecomplex *a, integer *lda, 
	integer *ipiv, doublecomplex *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zgetrs)(char *trans, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, integer *ipiv, doublecomplex *b, 
	integer *ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zggbak)(char *job, char *side, integer *n, integer *ilo, 
	integer *ihi, doublereal *lscale, doublereal *rscale, integer *m, 
	doublecomplex *v, integer *ldv, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zggbal)(char *job, integer *n, doublecomplex *a, integer 
	*lda, doublecomplex *b, integer *ldb, integer *ilo, integer *ihi, 
	doublereal *lscale, doublereal *rscale, doublereal *work, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(zgges)(char *jobvsl, char *jobvsr, char *sort, L_fp 
	delctg, integer *n, doublecomplex *a, integer *lda, doublecomplex *b, 
	integer *ldb, integer *sdim, doublecomplex *alpha, doublecomplex *
	beta, doublecomplex *vsl, integer *ldvsl, doublecomplex *vsr, integer 
	*ldvsr, doublecomplex *work, integer *lwork, doublereal *rwork, 
	logical *bwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zggesx)(char *jobvsl, char *jobvsr, char *sort, L_fp 
	delctg, char *sense, integer *n, doublecomplex *a, integer *lda, 
	doublecomplex *b, integer *ldb, integer *sdim, doublecomplex *alpha, 
	doublecomplex *beta, doublecomplex *vsl, integer *ldvsl, 
	doublecomplex *vsr, integer *ldvsr, doublereal *rconde, doublereal *
	rcondv, doublecomplex *work, integer *lwork, doublereal *rwork, 
	integer *iwork, integer *liwork, logical *bwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zggev)(char *jobvl, char *jobvr, integer *n, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublecomplex *alpha, doublecomplex *beta, doublecomplex *vl, integer 
	*ldvl, doublecomplex *vr, integer *ldvr, doublecomplex *work, integer 
	*lwork, doublereal *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zggevx)(char *balanc, char *jobvl, char *jobvr, char *
	sense, integer *n, doublecomplex *a, integer *lda, doublecomplex *b, 
	integer *ldb, doublecomplex *alpha, doublecomplex *beta, 
	doublecomplex *vl, integer *ldvl, doublecomplex *vr, integer *ldvr, 
	integer *ilo, integer *ihi, doublereal *lscale, doublereal *rscale, 
	doublereal *abnrm, doublereal *bbnrm, doublereal *rconde, doublereal *
	rcondv, doublecomplex *work, integer *lwork, doublereal *rwork, 
	integer *iwork, logical *bwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zggglm)(integer *n, integer *m, integer *p, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublecomplex *d__, doublecomplex *x, doublecomplex *y, doublecomplex 
	*work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zgghrd)(char *compq, char *compz, integer *n, integer *
	ilo, integer *ihi, doublecomplex *a, integer *lda, doublecomplex *b, 
	integer *ldb, doublecomplex *q, integer *ldq, doublecomplex *z__, 
	integer *ldz, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zgglse)(integer *m, integer *n, integer *p, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublecomplex *c__, doublecomplex *d__, doublecomplex *x, 
	doublecomplex *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zggqrf)(integer *n, integer *m, integer *p, 
	doublecomplex *a, integer *lda, doublecomplex *taua, doublecomplex *b,
	 integer *ldb, doublecomplex *taub, doublecomplex *work, integer *
	lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zggrqf)(integer *m, integer *p, integer *n, 
	doublecomplex *a, integer *lda, doublecomplex *taua, doublecomplex *b,
	 integer *ldb, doublecomplex *taub, doublecomplex *work, integer *
	lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zggsvd)(char *jobu, char *jobv, char *jobq, integer *m, 
	integer *n, integer *p, integer *k, integer *l, doublecomplex *a, 
	integer *lda, doublecomplex *b, integer *ldb, doublereal *alpha, 
	doublereal *beta, doublecomplex *u, integer *ldu, doublecomplex *v, 
	integer *ldv, doublecomplex *q, integer *ldq, doublecomplex *work, 
	doublereal *rwork, integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zggsvp)(char *jobu, char *jobv, char *jobq, integer *m, 
	integer *p, integer *n, doublecomplex *a, integer *lda, doublecomplex 
	*b, integer *ldb, doublereal *tola, doublereal *tolb, integer *k, 
	integer *l, doublecomplex *u, integer *ldu, doublecomplex *v, integer 
	*ldv, doublecomplex *q, integer *ldq, integer *iwork, doublereal *
	rwork, doublecomplex *tau, doublecomplex *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zgtcon)(char *norm, integer *n, doublecomplex *dl, 
	doublecomplex *d__, doublecomplex *du, doublecomplex *du2, integer *
	ipiv, doublereal *anorm, doublereal *rcond, doublecomplex *work, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zgtrfs)(char *trans, integer *n, integer *nrhs, 
	doublecomplex *dl, doublecomplex *d__, doublecomplex *du, 
	doublecomplex *dlf, doublecomplex *df, doublecomplex *duf, 
	doublecomplex *du2, integer *ipiv, doublecomplex *b, integer *ldb, 
	doublecomplex *x, integer *ldx, doublereal *ferr, doublereal *berr, 
	doublecomplex *work, doublereal *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zgtsv)(integer *n, integer *nrhs, doublecomplex *dl, 
	doublecomplex *d__, doublecomplex *du, doublecomplex *b, integer *ldb,
	 integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zgtsvx)(char *fact, char *trans, integer *n, integer *
	nrhs, doublecomplex *dl, doublecomplex *d__, doublecomplex *du, 
	doublecomplex *dlf, doublecomplex *df, doublecomplex *duf, 
	doublecomplex *du2, integer *ipiv, doublecomplex *b, integer *ldb, 
	doublecomplex *x, integer *ldx, doublereal *rcond, doublereal *ferr, 
	doublereal *berr, doublecomplex *work, doublereal *rwork, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(zgttrf)(integer *n, doublecomplex *dl, doublecomplex *
	d__, doublecomplex *du, doublecomplex *du2, integer *ipiv, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(zgttrs)(char *trans, integer *n, integer *nrhs, 
	doublecomplex *dl, doublecomplex *d__, doublecomplex *du, 
	doublecomplex *du2, integer *ipiv, doublecomplex *b, integer *ldb, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zgtts2)(integer *itrans, integer *n, integer *nrhs, 
	doublecomplex *dl, doublecomplex *d__, doublecomplex *du, 
	doublecomplex *du2, integer *ipiv, doublecomplex *b, integer *ldb);
 
/* Subroutine */ int LAPACKMANGLE(zhbev)(char *jobz, char *uplo, integer *n, integer *kd, 
	doublecomplex *ab, integer *ldab, doublereal *w, doublecomplex *z__, 
	integer *ldz, doublecomplex *work, doublereal *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zhbevd)(char *jobz, char *uplo, integer *n, integer *kd, 
	doublecomplex *ab, integer *ldab, doublereal *w, doublecomplex *z__, 
	integer *ldz, doublecomplex *work, integer *lwork, doublereal *rwork, 
	integer *lrwork, integer *iwork, integer *liwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zhbevx)(char *jobz, char *range, char *uplo, integer *n, 
	integer *kd, doublecomplex *ab, integer *ldab, doublecomplex *q, 
	integer *ldq, doublereal *vl, doublereal *vu, integer *il, integer *
	iu, doublereal *abstol, integer *m, doublereal *w, doublecomplex *z__,
	 integer *ldz, doublecomplex *work, doublereal *rwork, integer *iwork,
	 integer *ifail, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zhbgst)(char *vect, char *uplo, integer *n, integer *ka, 
	integer *kb, doublecomplex *ab, integer *ldab, doublecomplex *bb, 
	integer *ldbb, doublecomplex *x, integer *ldx, doublecomplex *work, 
	doublereal *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zhbgv)(char *jobz, char *uplo, integer *n, integer *ka, 
	integer *kb, doublecomplex *ab, integer *ldab, doublecomplex *bb, 
	integer *ldbb, doublereal *w, doublecomplex *z__, integer *ldz, 
	doublecomplex *work, doublereal *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zhbgvx)(char *jobz, char *range, char *uplo, integer *n, 
	integer *ka, integer *kb, doublecomplex *ab, integer *ldab, 
	doublecomplex *bb, integer *ldbb, doublecomplex *q, integer *ldq, 
	doublereal *vl, doublereal *vu, integer *il, integer *iu, doublereal *
	abstol, integer *m, doublereal *w, doublecomplex *z__, integer *ldz, 
	doublecomplex *work, doublereal *rwork, integer *iwork, integer *
	ifail, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zhbtrd)(char *vect, char *uplo, integer *n, integer *kd, 
	doublecomplex *ab, integer *ldab, doublereal *d__, doublereal *e, 
	doublecomplex *q, integer *ldq, doublecomplex *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zhecon)(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *ipiv, doublereal *anorm, doublereal *rcond, 
	doublecomplex *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zheev)(char *jobz, char *uplo, integer *n, doublecomplex 
	*a, integer *lda, doublereal *w, doublecomplex *work, integer *lwork, 
	doublereal *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zheevd)(char *jobz, char *uplo, integer *n, 
	doublecomplex *a, integer *lda, doublereal *w, doublecomplex *work, 
	integer *lwork, doublereal *rwork, integer *lrwork, integer *iwork, 
	integer *liwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zheevr)(char *jobz, char *range, char *uplo, integer *n, 
	doublecomplex *a, integer *lda, doublereal *vl, doublereal *vu, 
	integer *il, integer *iu, doublereal *abstol, integer *m, doublereal *
	w, doublecomplex *z__, integer *ldz, integer *isuppz, doublecomplex *
	work, integer *lwork, doublereal *rwork, integer *lrwork, integer *
	iwork, integer *liwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zheevx)(char *jobz, char *range, char *uplo, integer *n, 
	doublecomplex *a, integer *lda, doublereal *vl, doublereal *vu, 
	integer *il, integer *iu, doublereal *abstol, integer *m, doublereal *
	w, doublecomplex *z__, integer *ldz, doublecomplex *work, integer *
	lwork, doublereal *rwork, integer *iwork, integer *ifail, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(zhegs2)(integer *itype, char *uplo, integer *n, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zhegst)(integer *itype, char *uplo, integer *n, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zhegv)(integer *itype, char *jobz, char *uplo, integer *
	n, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublereal *w, doublecomplex *work, integer *lwork, doublereal *rwork,
	 integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zhegvd)(integer *itype, char *jobz, char *uplo, integer *
	n, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublereal *w, doublecomplex *work, integer *lwork, doublereal *rwork,
	 integer *lrwork, integer *iwork, integer *liwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zhegvx)(integer *itype, char *jobz, char *range, char *
	uplo, integer *n, doublecomplex *a, integer *lda, doublecomplex *b, 
	integer *ldb, doublereal *vl, doublereal *vu, integer *il, integer *
	iu, doublereal *abstol, integer *m, doublereal *w, doublecomplex *z__,
	 integer *ldz, doublecomplex *work, integer *lwork, doublereal *rwork,
	 integer *iwork, integer *ifail, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zherfs)(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, doublecomplex *af, integer *ldaf, 
	integer *ipiv, doublecomplex *b, integer *ldb, doublecomplex *x, 
	integer *ldx, doublereal *ferr, doublereal *berr, doublecomplex *work,
	 doublereal *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zhesv)(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, integer *ipiv, doublecomplex *b, 
	integer *ldb, doublecomplex *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zhesvx)(char *fact, char *uplo, integer *n, integer *
	nrhs, doublecomplex *a, integer *lda, doublecomplex *af, integer *
	ldaf, integer *ipiv, doublecomplex *b, integer *ldb, doublecomplex *x,
	 integer *ldx, doublereal *rcond, doublereal *ferr, doublereal *berr, 
	doublecomplex *work, integer *lwork, doublereal *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zhetf2)(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *ipiv, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zhetrd)(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, doublereal *d__, doublereal *e, doublecomplex *tau, 
	doublecomplex *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zhetrf)(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *ipiv, doublecomplex *work, integer *lwork, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zhetri)(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *ipiv, doublecomplex *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zhetrs)(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, integer *ipiv, doublecomplex *b, 
	integer *ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zhgeqz)(char *job, char *compq, char *compz, integer *n, 
	integer *ilo, integer *ihi, doublecomplex *a, integer *lda, 
	doublecomplex *b, integer *ldb, doublecomplex *alpha, doublecomplex *
	beta, doublecomplex *q, integer *ldq, doublecomplex *z__, integer *
	ldz, doublecomplex *work, integer *lwork, doublereal *rwork, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(zhpcon)(char *uplo, integer *n, doublecomplex *ap, 
	integer *ipiv, doublereal *anorm, doublereal *rcond, doublecomplex *
	work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zhpev)(char *jobz, char *uplo, integer *n, doublecomplex 
	*ap, doublereal *w, doublecomplex *z__, integer *ldz, doublecomplex *
	work, doublereal *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zhpevd)(char *jobz, char *uplo, integer *n, 
	doublecomplex *ap, doublereal *w, doublecomplex *z__, integer *ldz, 
	doublecomplex *work, integer *lwork, doublereal *rwork, integer *
	lrwork, integer *iwork, integer *liwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zhpevx)(char *jobz, char *range, char *uplo, integer *n, 
	doublecomplex *ap, doublereal *vl, doublereal *vu, integer *il, 
	integer *iu, doublereal *abstol, integer *m, doublereal *w, 
	doublecomplex *z__, integer *ldz, doublecomplex *work, doublereal *
	rwork, integer *iwork, integer *ifail, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zhpgst)(integer *itype, char *uplo, integer *n, 
	doublecomplex *ap, doublecomplex *bp, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zhpgv)(integer *itype, char *jobz, char *uplo, integer *
	n, doublecomplex *ap, doublecomplex *bp, doublereal *w, doublecomplex 
	*z__, integer *ldz, doublecomplex *work, doublereal *rwork, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(zhpgvd)(integer *itype, char *jobz, char *uplo, integer *
	n, doublecomplex *ap, doublecomplex *bp, doublereal *w, doublecomplex 
	*z__, integer *ldz, doublecomplex *work, integer *lwork, doublereal *
	rwork, integer *lrwork, integer *iwork, integer *liwork, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(zhpgvx)(integer *itype, char *jobz, char *range, char *
	uplo, integer *n, doublecomplex *ap, doublecomplex *bp, doublereal *
	vl, doublereal *vu, integer *il, integer *iu, doublereal *abstol, 
	integer *m, doublereal *w, doublecomplex *z__, integer *ldz, 
	doublecomplex *work, doublereal *rwork, integer *iwork, integer *
	ifail, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zhprfs)(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *ap, doublecomplex *afp, integer *ipiv, doublecomplex *
	b, integer *ldb, doublecomplex *x, integer *ldx, doublereal *ferr, 
	doublereal *berr, doublecomplex *work, doublereal *rwork, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(zhpsv)(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *ap, integer *ipiv, doublecomplex *b, integer *ldb, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zhpsvx)(char *fact, char *uplo, integer *n, integer *
	nrhs, doublecomplex *ap, doublecomplex *afp, integer *ipiv, 
	doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, 
	doublereal *rcond, doublereal *ferr, doublereal *berr, doublecomplex *
	work, doublereal *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zhptrd)(char *uplo, integer *n, doublecomplex *ap, 
	doublereal *d__, doublereal *e, doublecomplex *tau, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zhptrf)(char *uplo, integer *n, doublecomplex *ap, 
	integer *ipiv, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zhptri)(char *uplo, integer *n, doublecomplex *ap, 
	integer *ipiv, doublecomplex *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zhptrs)(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *ap, integer *ipiv, doublecomplex *b, integer *ldb, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zhsein)(char *side, char *eigsrc, char *initv, logical *
	select, integer *n, doublecomplex *h__, integer *ldh, doublecomplex *
	w, doublecomplex *vl, integer *ldvl, doublecomplex *vr, integer *ldvr,
	 integer *mm, integer *m, doublecomplex *work, doublereal *rwork, 
	integer *ifaill, integer *ifailr, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zhseqr)(char *job, char *compz, integer *n, integer *ilo,
	 integer *ihi, doublecomplex *h__, integer *ldh, doublecomplex *w, 
	doublecomplex *z__, integer *ldz, doublecomplex *work, integer *lwork,
	 integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zlabrd)(integer *m, integer *n, integer *nb, 
	doublecomplex *a, integer *lda, doublereal *d__, doublereal *e, 
	doublecomplex *tauq, doublecomplex *taup, doublecomplex *x, integer *
	ldx, doublecomplex *y, integer *ldy);
 
/* Subroutine */ int LAPACKMANGLE(zlacgv)(integer *n, doublecomplex *x, integer *incx);
 
/* Subroutine */ int LAPACKMANGLE(zlacon)(integer *n, doublecomplex *v, doublecomplex *x, 
	doublereal *est, integer *kase);
 
/* Subroutine */ int LAPACKMANGLE(zlacp2)(char *uplo, integer *m, integer *n, doublereal *
	a, integer *lda, doublecomplex *b, integer *ldb);
 
/* Subroutine */ int LAPACKMANGLE(zlacpy)(char *uplo, integer *m, integer *n, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb);
 
/* Subroutine */ int LAPACKMANGLE(zlacrm)(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublereal *b, integer *ldb, doublecomplex *c__, 
	integer *ldc, doublereal *rwork);
 
/* Subroutine */ int LAPACKMANGLE(zlacrt)(integer *n, doublecomplex *cx, integer *incx, 
	doublecomplex *cy, integer *incy, doublecomplex *c__, doublecomplex *
	s);
 
/* Subroutine */ int LAPACKMANGLE(zlaed0)(integer *qsiz, integer *n, doublereal *d__, 
	doublereal *e, doublecomplex *q, integer *ldq, doublecomplex *qstore, 
	integer *ldqs, doublereal *rwork, integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zlaed7)(integer *n, integer *cutpnt, integer *qsiz, 
	integer *tlvls, integer *curlvl, integer *curpbm, doublereal *d__, 
	doublecomplex *q, integer *ldq, doublereal *rho, integer *indxq, 
	doublereal *qstore, integer *qptr, integer *prmptr, integer *perm, 
	integer *givptr, integer *givcol, doublereal *givnum, doublecomplex *
	work, doublereal *rwork, integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zlaed8)(integer *k, integer *n, integer *qsiz, 
	doublecomplex *q, integer *ldq, doublereal *d__, doublereal *rho, 
	integer *cutpnt, doublereal *z__, doublereal *dlamda, doublecomplex *
	q2, integer *ldq2, doublereal *w, integer *indxp, integer *indx, 
	integer *indxq, integer *perm, integer *givptr, integer *givcol, 
	doublereal *givnum, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zlaein)(logical *rightv, logical *noinit, integer *n, 
	doublecomplex *h__, integer *ldh, doublecomplex *w, doublecomplex *v, 
	doublecomplex *b, integer *ldb, doublereal *rwork, doublereal *eps3, 
	doublereal *smlnum, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zlaesy)(doublecomplex *a, doublecomplex *b, 
	doublecomplex *c__, doublecomplex *rt1, doublecomplex *rt2, 
	doublecomplex *evscal, doublecomplex *cs1, doublecomplex *sn1);
 
/* Subroutine */ int LAPACKMANGLE(zlaev2)(doublecomplex *a, doublecomplex *b, 
	doublecomplex *c__, doublereal *rt1, doublereal *rt2, doublereal *cs1,
	 doublecomplex *sn1);
 
/* Subroutine */ int LAPACKMANGLE(zlags2)(logical *upper, doublereal *a1, doublecomplex *
	a2, doublereal *a3, doublereal *b1, doublecomplex *b2, doublereal *b3,
	 doublereal *csu, doublecomplex *snu, doublereal *csv, doublecomplex *
	snv, doublereal *csq, doublecomplex *snq);
 
/* Subroutine */ int LAPACKMANGLE(zlagtm)(char *trans, integer *n, integer *nrhs, 
	doublereal *alpha, doublecomplex *dl, doublecomplex *d__, 
	doublecomplex *du, doublecomplex *x, integer *ldx, doublereal *beta, 
	doublecomplex *b, integer *ldb);
 
/* Subroutine */ int LAPACKMANGLE(zlahef)(char *uplo, integer *n, integer *nb, integer *kb,
	 doublecomplex *a, integer *lda, integer *ipiv, doublecomplex *w, 
	integer *ldw, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zlahqr)(logical *wantt, logical *wantz, integer *n, 
	integer *ilo, integer *ihi, doublecomplex *h__, integer *ldh, 
	doublecomplex *w, integer *iloz, integer *ihiz, doublecomplex *z__, 
	integer *ldz, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zlahrd)(integer *n, integer *k, integer *nb, 
	doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *t, 
	integer *ldt, doublecomplex *y, integer *ldy);
 
/* Subroutine */ int LAPACKMANGLE(zlaic1)(integer *job, integer *j, doublecomplex *x, 
	doublereal *sest, doublecomplex *w, doublecomplex *gamma, doublereal *
	sestpr, doublecomplex *s, doublecomplex *c__);
 
/* Subroutine */ int LAPACKMANGLE(zlals0)(integer *icompq, integer *nl, integer *nr, 
	integer *sqre, integer *nrhs, doublecomplex *b, integer *ldb, 
	doublecomplex *bx, integer *ldbx, integer *perm, integer *givptr, 
	integer *givcol, integer *ldgcol, doublereal *givnum, integer *ldgnum,
	 doublereal *poles, doublereal *difl, doublereal *difr, doublereal *
	z__, integer *k, doublereal *c__, doublereal *s, doublereal *rwork, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zlalsa)(integer *icompq, integer *smlsiz, integer *n, 
	integer *nrhs, doublecomplex *b, integer *ldb, doublecomplex *bx, 
	integer *ldbx, doublereal *u, integer *ldu, doublereal *vt, integer *
	k, doublereal *difl, doublereal *difr, doublereal *z__, doublereal *
	poles, integer *givptr, integer *givcol, integer *ldgcol, integer *
	perm, doublereal *givnum, doublereal *c__, doublereal *s, doublereal *
	rwork, integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zlapll)(integer *n, doublecomplex *x, integer *incx, 
	doublecomplex *y, integer *incy, doublereal *ssmin);
 
/* Subroutine */ int LAPACKMANGLE(zlapmt)(logical *forwrd, integer *m, integer *n, 
	doublecomplex *x, integer *ldx, integer *k);
 
/* Subroutine */ int LAPACKMANGLE(zlaqgb)(integer *m, integer *n, integer *kl, integer *ku,
	 doublecomplex *ab, integer *ldab, doublereal *r__, doublereal *c__, 
	doublereal *rowcnd, doublereal *colcnd, doublereal *amax, char *equed);
 
/* Subroutine */ int LAPACKMANGLE(zlaqge)(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublereal *r__, doublereal *c__, doublereal *rowcnd, 
	doublereal *colcnd, doublereal *amax, char *equed);
 
/* Subroutine */ int LAPACKMANGLE(zlaqhb)(char *uplo, integer *n, integer *kd, 
	doublecomplex *ab, integer *ldab, doublereal *s, doublereal *scond, 
	doublereal *amax, char *equed);
 
/* Subroutine */ int LAPACKMANGLE(zlaqhe)(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, doublereal *s, doublereal *scond, doublereal *amax, 
	char *equed);
 
/* Subroutine */ int LAPACKMANGLE(zlaqhp)(char *uplo, integer *n, doublecomplex *ap, 
	doublereal *s, doublereal *scond, doublereal *amax, char *equed);
 
/* Subroutine */ int LAPACKMANGLE(zlaqp2)(integer *m, integer *n, integer *offset, 
	doublecomplex *a, integer *lda, integer *jpvt, doublecomplex *tau, 
	doublereal *vn1, doublereal *vn2, doublecomplex *work);
 
/* Subroutine */ int LAPACKMANGLE(zlaqps)(integer *m, integer *n, integer *offset, integer 
	*nb, integer *kb, doublecomplex *a, integer *lda, integer *jpvt, 
	doublecomplex *tau, doublereal *vn1, doublereal *vn2, doublecomplex *
	auxv, doublecomplex *f, integer *ldf);
 
/* Subroutine */ int LAPACKMANGLE(zlaqsb)(char *uplo, integer *n, integer *kd, 
	doublecomplex *ab, integer *ldab, doublereal *s, doublereal *scond, 
	doublereal *amax, char *equed);
 
/* Subroutine */ int LAPACKMANGLE(zlaqsp)(char *uplo, integer *n, doublecomplex *ap, 
	doublereal *s, doublereal *scond, doublereal *amax, char *equed);
 
/* Subroutine */ int LAPACKMANGLE(zlaqsy)(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, doublereal *s, doublereal *scond, doublereal *amax, 
	char *equed);
 
/* Subroutine */ int LAPACKMANGLE(zlar1v)(integer *n, integer *b1, integer *bn, doublereal 
	*sigma, doublereal *d__, doublereal *l, doublereal *ld, doublereal *
	lld, doublereal *gersch, doublecomplex *z__, doublereal *ztz, 
	doublereal *mingma, integer *r__, integer *isuppz, doublereal *work);
 
/* Subroutine */ int LAPACKMANGLE(zlar2v)(integer *n, doublecomplex *x, doublecomplex *y, 
	doublecomplex *z__, integer *incx, doublereal *c__, doublecomplex *s, 
	integer *incc);
 
/* Subroutine */ int LAPACKMANGLE(zlarcm)(integer *m, integer *n, doublereal *a, integer *
	lda, doublecomplex *b, integer *ldb, doublecomplex *c__, integer *ldc,
	 doublereal *rwork);
 
/* Subroutine */ int LAPACKMANGLE(zlarf)(char *side, integer *m, integer *n, doublecomplex 
	*v, integer *incv, doublecomplex *tau, doublecomplex *c__, integer *
	ldc, doublecomplex *work);
 
/* Subroutine */ int LAPACKMANGLE(zlarfb)(char *side, char *trans, char *direct, char *
	storev, integer *m, integer *n, integer *k, doublecomplex *v, integer 
	*ldv, doublecomplex *t, integer *ldt, doublecomplex *c__, integer *
	ldc, doublecomplex *work, integer *ldwork);
 
/* Subroutine */ int LAPACKMANGLE(zlarfg)(integer *n, doublecomplex *alpha, doublecomplex *
	x, integer *incx, doublecomplex *tau);
 
/* Subroutine */ int LAPACKMANGLE(zlarft)(char *direct, char *storev, integer *n, integer *
	k, doublecomplex *v, integer *ldv, doublecomplex *tau, doublecomplex *
	t, integer *ldt);
 
/* Subroutine */ int LAPACKMANGLE(zlarfx)(char *side, integer *m, integer *n, 
	doublecomplex *v, doublecomplex *tau, doublecomplex *c__, integer *
	ldc, doublecomplex *work);
 
/* Subroutine */ int LAPACKMANGLE(zlargv)(integer *n, doublecomplex *x, integer *incx, 
	doublecomplex *y, integer *incy, doublereal *c__, integer *incc);
 
/* Subroutine */ int LAPACKMANGLE(zlarnv)(integer *idist, integer *iseed, integer *n, 
	doublecomplex *x);
 
/* Subroutine */ int LAPACKMANGLE(zlarrv)(integer *n, doublereal *d__, doublereal *l, 
	integer *isplit, integer *m, doublereal *w, integer *iblock, 
	doublereal *gersch, doublereal *tol, doublecomplex *z__, integer *ldz,
	 integer *isuppz, doublereal *work, integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zlartg)(doublecomplex *f, doublecomplex *g, doublereal *
	cs, doublecomplex *sn, doublecomplex *r__);
 
/* Subroutine */ int LAPACKMANGLE(zlartv)(integer *n, doublecomplex *x, integer *incx, 
	doublecomplex *y, integer *incy, doublereal *c__, doublecomplex *s, 
	integer *incc);
 
/* Subroutine */ int LAPACKMANGLE(zlarz)(char *side, integer *m, integer *n, integer *l, 
	doublecomplex *v, integer *incv, doublecomplex *tau, doublecomplex *
	c__, integer *ldc, doublecomplex *work);
 
/* Subroutine */ int LAPACKMANGLE(zlarzb)(char *side, char *trans, char *direct, char *
	storev, integer *m, integer *n, integer *k, integer *l, doublecomplex 
	*v, integer *ldv, doublecomplex *t, integer *ldt, doublecomplex *c__, 
	integer *ldc, doublecomplex *work, integer *ldwork);
 
/* Subroutine */ int LAPACKMANGLE(zlarzt)(char *direct, char *storev, integer *n, integer *
	k, doublecomplex *v, integer *ldv, doublecomplex *tau, doublecomplex *
	t, integer *ldt);
 
/* Subroutine */ int LAPACKMANGLE(zlascl)(char *type__, integer *kl, integer *ku, 
	doublereal *cfrom, doublereal *cto, integer *m, integer *n, 
	doublecomplex *a, integer *lda, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zlaset)(char *uplo, integer *m, integer *n, 
	doublecomplex *alpha, doublecomplex *beta, doublecomplex *a, integer *
	lda);
 
/* Subroutine */ int LAPACKMANGLE(zlasr)(char *side, char *pivot, char *direct, integer *m,
	 integer *n, doublereal *c__, doublereal *s, doublecomplex *a, 
	integer *lda);
 
/* Subroutine */ int LAPACKMANGLE(zlassq)(integer *n, doublecomplex *x, integer *incx, 
	doublereal *scale, doublereal *sumsq);
 
/* Subroutine */ int LAPACKMANGLE(zlaswp)(integer *n, doublecomplex *a, integer *lda, 
	integer *k1, integer *k2, integer *ipiv, integer *incx);
 
/* Subroutine */ int LAPACKMANGLE(zlasyf)(char *uplo, integer *n, integer *nb, integer *kb,
	 doublecomplex *a, integer *lda, integer *ipiv, doublecomplex *w, 
	integer *ldw, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zlatbs)(char *uplo, char *trans, char *diag, char *
	normin, integer *n, integer *kd, doublecomplex *ab, integer *ldab, 
	doublecomplex *x, doublereal *scale, doublereal *cnorm, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zlatdf)(integer *ijob, integer *n, doublecomplex *z__, 
	integer *ldz, doublecomplex *rhs, doublereal *rdsum, doublereal *
	rdscal, integer *ipiv, integer *jpiv);
 
/* Subroutine */ int LAPACKMANGLE(zlatps)(char *uplo, char *trans, char *diag, char *
	normin, integer *n, doublecomplex *ap, doublecomplex *x, doublereal *
	scale, doublereal *cnorm, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zlatrd)(char *uplo, integer *n, integer *nb, 
	doublecomplex *a, integer *lda, doublereal *e, doublecomplex *tau, 
	doublecomplex *w, integer *ldw);
 
/* Subroutine */ int LAPACKMANGLE(zlatrs)(char *uplo, char *trans, char *diag, char *
	normin, integer *n, doublecomplex *a, integer *lda, doublecomplex *x, 
	doublereal *scale, doublereal *cnorm, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zlatrz)(integer *m, integer *n, integer *l, 
	doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
	work);
 
/* Subroutine */ int LAPACKMANGLE(zlatzm)(char *side, integer *m, integer *n, 
	doublecomplex *v, integer *incv, doublecomplex *tau, doublecomplex *
	c1, doublecomplex *c2, integer *ldc, doublecomplex *work);
 
/* Subroutine */ int LAPACKMANGLE(zlauu2)(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zlauum)(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zpbcon)(char *uplo, integer *n, integer *kd, 
	doublecomplex *ab, integer *ldab, doublereal *anorm, doublereal *
	rcond, doublecomplex *work, doublereal *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zpbequ)(char *uplo, integer *n, integer *kd, 
	doublecomplex *ab, integer *ldab, doublereal *s, doublereal *scond, 
	doublereal *amax, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zpbrfs)(char *uplo, integer *n, integer *kd, integer *
	nrhs, doublecomplex *ab, integer *ldab, doublecomplex *afb, integer *
	ldafb, doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx,
	 doublereal *ferr, doublereal *berr, doublecomplex *work, doublereal *
	rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zpbstf)(char *uplo, integer *n, integer *kd, 
	doublecomplex *ab, integer *ldab, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zpbsv)(char *uplo, integer *n, integer *kd, integer *
	nrhs, doublecomplex *ab, integer *ldab, doublecomplex *b, integer *
	ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zpbsvx)(char *fact, char *uplo, integer *n, integer *kd, 
	integer *nrhs, doublecomplex *ab, integer *ldab, doublecomplex *afb, 
	integer *ldafb, char *equed, doublereal *s, doublecomplex *b, integer 
	*ldb, doublecomplex *x, integer *ldx, doublereal *rcond, doublereal *
	ferr, doublereal *berr, doublecomplex *work, doublereal *rwork, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zpbtf2)(char *uplo, integer *n, integer *kd, 
	doublecomplex *ab, integer *ldab, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zpbtrf)(char *uplo, integer *n, integer *kd, 
	doublecomplex *ab, integer *ldab, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zpbtrs)(char *uplo, integer *n, integer *kd, integer *
	nrhs, doublecomplex *ab, integer *ldab, doublecomplex *b, integer *
	ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zpocon)(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, doublereal *anorm, doublereal *rcond, doublecomplex *
	work, doublereal *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zpoequ)(integer *n, doublecomplex *a, integer *lda, 
	doublereal *s, doublereal *scond, doublereal *amax, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zporfs)(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, doublecomplex *af, integer *ldaf, 
	doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, 
	doublereal *ferr, doublereal *berr, doublecomplex *work, doublereal *
	rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zposv)(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zposvx)(char *fact, char *uplo, integer *n, integer *
	nrhs, doublecomplex *a, integer *lda, doublecomplex *af, integer *
	ldaf, char *equed, doublereal *s, doublecomplex *b, integer *ldb, 
	doublecomplex *x, integer *ldx, doublereal *rcond, doublereal *ferr, 
	doublereal *berr, doublecomplex *work, doublereal *rwork, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(zpotf2)(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zpotrf)(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zpotri)(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zpotrs)(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zppcon)(char *uplo, integer *n, doublecomplex *ap, 
	doublereal *anorm, doublereal *rcond, doublecomplex *work, doublereal 
	*rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zppequ)(char *uplo, integer *n, doublecomplex *ap, 
	doublereal *s, doublereal *scond, doublereal *amax, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zpprfs)(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *ap, doublecomplex *afp, doublecomplex *b, integer *ldb,
	 doublecomplex *x, integer *ldx, doublereal *ferr, doublereal *berr, 
	doublecomplex *work, doublereal *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zppsv)(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *ap, doublecomplex *b, integer *ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zppsvx)(char *fact, char *uplo, integer *n, integer *
	nrhs, doublecomplex *ap, doublecomplex *afp, char *equed, doublereal *
	s, doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, 
	doublereal *rcond, doublereal *ferr, doublereal *berr, doublecomplex *
	work, doublereal *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zpptrf)(char *uplo, integer *n, doublecomplex *ap, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zpptri)(char *uplo, integer *n, doublecomplex *ap, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zpptrs)(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *ap, doublecomplex *b, integer *ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zptcon)(integer *n, doublereal *d__, doublecomplex *e, 
	doublereal *anorm, doublereal *rcond, doublereal *rwork, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(zptrfs)(char *uplo, integer *n, integer *nrhs, 
	doublereal *d__, doublecomplex *e, doublereal *df, doublecomplex *ef, 
	doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, 
	doublereal *ferr, doublereal *berr, doublecomplex *work, doublereal *
	rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zptsv)(integer *n, integer *nrhs, doublereal *d__, 
	doublecomplex *e, doublecomplex *b, integer *ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zptsvx)(char *fact, integer *n, integer *nrhs, 
	doublereal *d__, doublecomplex *e, doublereal *df, doublecomplex *ef, 
	doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, 
	doublereal *rcond, doublereal *ferr, doublereal *berr, doublecomplex *
	work, doublereal *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zpttrf)(integer *n, doublereal *d__, doublecomplex *e, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zpttrs)(char *uplo, integer *n, integer *nrhs, 
	doublereal *d__, doublecomplex *e, doublecomplex *b, integer *ldb, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zptts2)(integer *iuplo, integer *n, integer *nrhs, 
	doublereal *d__, doublecomplex *e, doublecomplex *b, integer *ldb);
 
/* Subroutine */ int LAPACKMANGLE(zrot)(integer *n, doublecomplex *cx, integer *incx, 
	doublecomplex *cy, integer *incy, doublereal *c__, doublecomplex *s);
 
/* Subroutine */ int LAPACKMANGLE(zspcon)(char *uplo, integer *n, doublecomplex *ap, 
	integer *ipiv, doublereal *anorm, doublereal *rcond, doublecomplex *
	work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zspmv)(char *uplo, integer *n, doublecomplex *alpha, 
	doublecomplex *ap, doublecomplex *x, integer *incx, doublecomplex *
	beta, doublecomplex *y, integer *incy);
 
/* Subroutine */ int LAPACKMANGLE(zspr)(char *uplo, integer *n, doublecomplex *alpha, 
	doublecomplex *x, integer *incx, doublecomplex *ap);
 
/* Subroutine */ int LAPACKMANGLE(zsprfs)(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *ap, doublecomplex *afp, integer *ipiv, doublecomplex *
	b, integer *ldb, doublecomplex *x, integer *ldx, doublereal *ferr, 
	doublereal *berr, doublecomplex *work, doublereal *rwork, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(zspsv)(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *ap, integer *ipiv, doublecomplex *b, integer *ldb, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zspsvx)(char *fact, char *uplo, integer *n, integer *
	nrhs, doublecomplex *ap, doublecomplex *afp, integer *ipiv, 
	doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, 
	doublereal *rcond, doublereal *ferr, doublereal *berr, doublecomplex *
	work, doublereal *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zsptrf)(char *uplo, integer *n, doublecomplex *ap, 
	integer *ipiv, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zsptri)(char *uplo, integer *n, doublecomplex *ap, 
	integer *ipiv, doublecomplex *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zsptrs)(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *ap, integer *ipiv, doublecomplex *b, integer *ldb, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zstedc)(char *compz, integer *n, doublereal *d__, 
	doublereal *e, doublecomplex *z__, integer *ldz, doublecomplex *work, 
	integer *lwork, doublereal *rwork, integer *lrwork, integer *iwork, 
	integer *liwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zstein)(integer *n, doublereal *d__, doublereal *e, 
	integer *m, doublereal *w, integer *iblock, integer *isplit, 
	doublecomplex *z__, integer *ldz, doublereal *work, integer *iwork, 
	integer *ifail, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zsteqr)(char *compz, integer *n, doublereal *d__, 
	doublereal *e, doublecomplex *z__, integer *ldz, doublereal *work, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zsycon)(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *ipiv, doublereal *anorm, doublereal *rcond, 
	doublecomplex *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zsymv)(char *uplo, integer *n, doublecomplex *alpha, 
	doublecomplex *a, integer *lda, doublecomplex *x, integer *incx, 
	doublecomplex *beta, doublecomplex *y, integer *incy);
 
/* Subroutine */ int LAPACKMANGLE(zsyr)(char *uplo, integer *n, doublecomplex *alpha, 
	doublecomplex *x, integer *incx, doublecomplex *a, integer *lda);
 
/* Subroutine */ int LAPACKMANGLE(zsyrfs)(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, doublecomplex *af, integer *ldaf, 
	integer *ipiv, doublecomplex *b, integer *ldb, doublecomplex *x, 
	integer *ldx, doublereal *ferr, doublereal *berr, doublecomplex *work,
	 doublereal *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zsysv)(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, integer *ipiv, doublecomplex *b, 
	integer *ldb, doublecomplex *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zsysvx)(char *fact, char *uplo, integer *n, integer *
	nrhs, doublecomplex *a, integer *lda, doublecomplex *af, integer *
	ldaf, integer *ipiv, doublecomplex *b, integer *ldb, doublecomplex *x,
	 integer *ldx, doublereal *rcond, doublereal *ferr, doublereal *berr, 
	doublecomplex *work, integer *lwork, doublereal *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zsytf2)(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *ipiv, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zsytrf)(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *ipiv, doublecomplex *work, integer *lwork, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zsytri)(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *ipiv, doublecomplex *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zsytrs)(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, integer *ipiv, doublecomplex *b, 
	integer *ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ztbcon)(char *norm, char *uplo, char *diag, integer *n, 
	integer *kd, doublecomplex *ab, integer *ldab, doublereal *rcond, 
	doublecomplex *work, doublereal *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ztbrfs)(char *uplo, char *trans, char *diag, integer *n, 
	integer *kd, integer *nrhs, doublecomplex *ab, integer *ldab, 
	doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, 
	doublereal *ferr, doublereal *berr, doublecomplex *work, doublereal *
	rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ztbtrs)(char *uplo, char *trans, char *diag, integer *n, 
	integer *kd, integer *nrhs, doublecomplex *ab, integer *ldab, 
	doublecomplex *b, integer *ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ztgevc)(char *side, char *howmny, logical *select, 
	integer *n, doublecomplex *a, integer *lda, doublecomplex *b, integer 
	*ldb, doublecomplex *vl, integer *ldvl, doublecomplex *vr, integer *
	ldvr, integer *mm, integer *m, doublecomplex *work, doublereal *rwork,
	 integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ztgex2)(logical *wantq, logical *wantz, integer *n, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublecomplex *q, integer *ldq, doublecomplex *z__, integer *ldz, 
	integer *j1, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ztgexc)(logical *wantq, logical *wantz, integer *n, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublecomplex *q, integer *ldq, doublecomplex *z__, integer *ldz, 
	integer *ifst, integer *ilst, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ztgsen)(integer *ijob, logical *wantq, logical *wantz, 
	logical *select, integer *n, doublecomplex *a, integer *lda, 
	doublecomplex *b, integer *ldb, doublecomplex *alpha, doublecomplex *
	beta, doublecomplex *q, integer *ldq, doublecomplex *z__, integer *
	ldz, integer *m, doublereal *pl, doublereal *pr, doublereal *dif, 
	doublecomplex *work, integer *lwork, integer *iwork, integer *liwork, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ztgsja)(char *jobu, char *jobv, char *jobq, integer *m, 
	integer *p, integer *n, integer *k, integer *l, doublecomplex *a, 
	integer *lda, doublecomplex *b, integer *ldb, doublereal *tola, 
	doublereal *tolb, doublereal *alpha, doublereal *beta, doublecomplex *
	u, integer *ldu, doublecomplex *v, integer *ldv, doublecomplex *q, 
	integer *ldq, doublecomplex *work, integer *ncycle, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ztgsna)(char *job, char *howmny, logical *select, 
	integer *n, doublecomplex *a, integer *lda, doublecomplex *b, integer 
	*ldb, doublecomplex *vl, integer *ldvl, doublecomplex *vr, integer *
	ldvr, doublereal *s, doublereal *dif, integer *mm, integer *m, 
	doublecomplex *work, integer *lwork, integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ztgsy2)(char *trans, integer *ijob, integer *m, integer *
	n, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublecomplex *c__, integer *ldc, doublecomplex *d__, integer *ldd, 
	doublecomplex *e, integer *lde, doublecomplex *f, integer *ldf, 
	doublereal *scale, doublereal *rdsum, doublereal *rdscal, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(ztgsyl)(char *trans, integer *ijob, integer *m, integer *
	n, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublecomplex *c__, integer *ldc, doublecomplex *d__, integer *ldd, 
	doublecomplex *e, integer *lde, doublecomplex *f, integer *ldf, 
	doublereal *scale, doublereal *dif, doublecomplex *work, integer *
	lwork, integer *iwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ztpcon)(char *norm, char *uplo, char *diag, integer *n, 
	doublecomplex *ap, doublereal *rcond, doublecomplex *work, doublereal 
	*rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ztprfs)(char *uplo, char *trans, char *diag, integer *n, 
	integer *nrhs, doublecomplex *ap, doublecomplex *b, integer *ldb, 
	doublecomplex *x, integer *ldx, doublereal *ferr, doublereal *berr, 
	doublecomplex *work, doublereal *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ztptri)(char *uplo, char *diag, integer *n, 
	doublecomplex *ap, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ztptrs)(char *uplo, char *trans, char *diag, integer *n, 
	integer *nrhs, doublecomplex *ap, doublecomplex *b, integer *ldb, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ztrcon)(char *norm, char *uplo, char *diag, integer *n, 
	doublecomplex *a, integer *lda, doublereal *rcond, doublecomplex *
	work, doublereal *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ztrevc)(char *side, char *howmny, logical *select, 
	integer *n, doublecomplex *t, integer *ldt, doublecomplex *vl, 
	integer *ldvl, doublecomplex *vr, integer *ldvr, integer *mm, integer 
	*m, doublecomplex *work, doublereal *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ztrexc)(char *compq, integer *n, doublecomplex *t, 
	integer *ldt, doublecomplex *q, integer *ldq, integer *ifst, integer *
	ilst, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ztrrfs)(char *uplo, char *trans, char *diag, integer *n, 
	integer *nrhs, doublecomplex *a, integer *lda, doublecomplex *b, 
	integer *ldb, doublecomplex *x, integer *ldx, doublereal *ferr, 
	doublereal *berr, doublecomplex *work, doublereal *rwork, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(ztrsen)(char *job, char *compq, logical *select, integer 
	*n, doublecomplex *t, integer *ldt, doublecomplex *q, integer *ldq, 
	doublecomplex *w, integer *m, doublereal *s, doublereal *sep, 
	doublecomplex *work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ztrsna)(char *job, char *howmny, logical *select, 
	integer *n, doublecomplex *t, integer *ldt, doublecomplex *vl, 
	integer *ldvl, doublecomplex *vr, integer *ldvr, doublereal *s, 
	doublereal *sep, integer *mm, integer *m, doublecomplex *work, 
	integer *ldwork, doublereal *rwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ztrsyl)(char *trana, char *tranb, integer *isgn, integer 
	*m, integer *n, doublecomplex *a, integer *lda, doublecomplex *b, 
	integer *ldb, doublecomplex *c__, integer *ldc, doublereal *scale, 
	integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ztrti2)(char *uplo, char *diag, integer *n, 
	doublecomplex *a, integer *lda, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ztrtri)(char *uplo, char *diag, integer *n, 
	doublecomplex *a, integer *lda, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ztrtrs)(char *uplo, char *trans, char *diag, integer *n, 
	integer *nrhs, doublecomplex *a, integer *lda, doublecomplex *b, 
	integer *ldb, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ztzrqf)(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublecomplex *tau, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(ztzrzf)(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublecomplex *tau, doublecomplex *work, integer *lwork,
	 integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zung2l)(integer *m, integer *n, integer *k, 
	doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
	work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zung2r)(integer *m, integer *n, integer *k, 
	doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
	work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zungbr)(char *vect, integer *m, integer *n, integer *k, 
	doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
	work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zunghr)(integer *n, integer *ilo, integer *ihi, 
	doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
	work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zungl2)(integer *m, integer *n, integer *k, 
	doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
	work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zunglq)(integer *m, integer *n, integer *k, 
	doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
	work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zungql)(integer *m, integer *n, integer *k, 
	doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
	work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zungqr)(integer *m, integer *n, integer *k, 
	doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
	work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zungr2)(integer *m, integer *n, integer *k, 
	doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
	work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zungrq)(integer *m, integer *n, integer *k, 
	doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
	work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zungtr)(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, doublecomplex *tau, doublecomplex *work, integer *lwork,
	 integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zunm2l)(char *side, char *trans, integer *m, integer *n, 
	integer *k, doublecomplex *a, integer *lda, doublecomplex *tau, 
	doublecomplex *c__, integer *ldc, doublecomplex *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zunm2r)(char *side, char *trans, integer *m, integer *n, 
	integer *k, doublecomplex *a, integer *lda, doublecomplex *tau, 
	doublecomplex *c__, integer *ldc, doublecomplex *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zunmbr)(char *vect, char *side, char *trans, integer *m, 
	integer *n, integer *k, doublecomplex *a, integer *lda, doublecomplex 
	*tau, doublecomplex *c__, integer *ldc, doublecomplex *work, integer *
	lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zunmhr)(char *side, char *trans, integer *m, integer *n, 
	integer *ilo, integer *ihi, doublecomplex *a, integer *lda, 
	doublecomplex *tau, doublecomplex *c__, integer *ldc, doublecomplex *
	work, integer *lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zunml2)(char *side, char *trans, integer *m, integer *n, 
	integer *k, doublecomplex *a, integer *lda, doublecomplex *tau, 
	doublecomplex *c__, integer *ldc, doublecomplex *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zunmlq)(char *side, char *trans, integer *m, integer *n, 
	integer *k, doublecomplex *a, integer *lda, doublecomplex *tau, 
	doublecomplex *c__, integer *ldc, doublecomplex *work, integer *lwork,
	 integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zunmql)(char *side, char *trans, integer *m, integer *n, 
	integer *k, doublecomplex *a, integer *lda, doublecomplex *tau, 
	doublecomplex *c__, integer *ldc, doublecomplex *work, integer *lwork,
	 integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zunmqr)(char *side, char *trans, integer *m, integer *n, 
	integer *k, doublecomplex *a, integer *lda, doublecomplex *tau, 
	doublecomplex *c__, integer *ldc, doublecomplex *work, integer *lwork,
	 integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zunmr2)(char *side, char *trans, integer *m, integer *n, 
	integer *k, doublecomplex *a, integer *lda, doublecomplex *tau, 
	doublecomplex *c__, integer *ldc, doublecomplex *work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zunmr3)(char *side, char *trans, integer *m, integer *n, 
	integer *k, integer *l, doublecomplex *a, integer *lda, doublecomplex 
	*tau, doublecomplex *c__, integer *ldc, doublecomplex *work, integer *
	info);
 
/* Subroutine */ int LAPACKMANGLE(zunmrq)(char *side, char *trans, integer *m, integer *n, 
	integer *k, doublecomplex *a, integer *lda, doublecomplex *tau, 
	doublecomplex *c__, integer *ldc, doublecomplex *work, integer *lwork,
	 integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zunmrz)(char *side, char *trans, integer *m, integer *n, 
	integer *k, integer *l, doublecomplex *a, integer *lda, doublecomplex 
	*tau, doublecomplex *c__, integer *ldc, doublecomplex *work, integer *
	lwork, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zunmtr)(char *side, char *uplo, char *trans, integer *m, 
	integer *n, doublecomplex *a, integer *lda, doublecomplex *tau, 
	doublecomplex *c__, integer *ldc, doublecomplex *work, integer *lwork,
	 integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zupgtr)(char *uplo, integer *n, doublecomplex *ap, 
	doublecomplex *tau, doublecomplex *q, integer *ldq, doublecomplex *
	work, integer *info);
 
/* Subroutine */ int LAPACKMANGLE(zupmtr)(char *side, char *uplo, char *trans, integer *m, 
	integer *n, doublecomplex *ap, doublecomplex *tau, doublecomplex *c__,
	 integer *ldc, doublecomplex *work, integer *info);
}
#endif /* __CLAPACK_H */
