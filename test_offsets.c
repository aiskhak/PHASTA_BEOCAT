// test_offsets.c
#include <stdio.h>
#include <stddef.h>
#include "common_c_structs.h"

#define OFF(f) printf("%-24s : %4zu\n", #f, offsetof(struct levlset, f))

int main(void) {
    OFF(epsilon_ls);
    OFF(epsilonBT);
    OFF(rholiq);
    OFF(coalbubrad);
	OFF(i_res_cf);
    OFF(iBT);
    OFF(coalcon_rem);
    printf("TOTAL SIZE               : %4zu\n", sizeof(struct levlset));
    return 0;
}
