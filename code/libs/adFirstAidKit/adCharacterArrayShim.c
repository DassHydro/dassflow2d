/*
 * adMPI.c calls pushCharacterArray()/popCharacterArray() directly as plain C
 * functions (mixed case, no trailing underscore -- see its #include of
 * adStack.h). DassFlow's existing push/pop stack (src/adjoint/adBuffer.f)
 * already implements the equivalent primitive, but only reachable through
 * gfortran's external-procedure name mangling (lowercase + trailing
 * underscore, passed by reference): pushcharacterarray_/popcharacterarray_.
 *
 * Rather than vendoring a second, conflicting stack implementation (see the
 * comment on df_sum_r in src/common/m_mpi.f90 for why a second adStack.c
 * cannot coexist with adBuffer.f), this is a thin bridge from the C-style
 * names adMPI.c needs onto the Fortran ones DassFlow already has and uses
 * everywhere else.
 */

extern void pushcharacterarray_(char *x, int *n);
extern void popcharacterarray_(char *x, int *n);
extern void pushinteger4_(int *x);
extern void popinteger4_(int *x);

void pushCharacterArray(char *x, int n) {
    pushcharacterarray_(x, &n);
}

void popCharacterArray(char *x, int n) {
    popcharacterarray_(x, &n);
}

void pushInteger4(int val) {
    pushinteger4_(&val);
}

void popInteger4(int *val) {
    popinteger4_(val);
}
