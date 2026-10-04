/* Count the terms of A367620 from a(1) = 20 to its first branch point 19999999918 one step at a time,
 * without block jumps, as an independent check of analysis/first_branch_point.py.
 *
 * Run with: gcc -O2 -o first_branch_point analysis/first_branch_point.c && ./first_branch_point
 */
#include <stdio.h>
#include <stdint.h>

static int leading(int64_t n) { while (n >= 10) n /= 10; return (int)n; }

int main(void) {
    const int64_t target = 19999999918LL;
    int64_t n = 20, index = 1, branch_points_before = 0;   /* a(1) = 20 */
    while (n != target) {
        int64_t x = n % 10, next = -1;
        int children = 0;
        for (int e = 1; e <= 9; e++) {
            int64_t c = n + 10 * x + e;
            if (leading(c) == e) { if (next < 0) next = c; children++; }
        }
        if (next < 0) { printf("sequence died at %lld (index %lld)\n", (long long)n, (long long)index); return 1; }
        if (children == 2) branch_points_before++;
        n = next;
        index++;
    }
    int64_t x = n % 10; int children = 0;
    for (int e = 1; e <= 9; e++) if (leading(n + 10 * x + e) == e) children++;
    printf("%lld is term %lld (a(1) = 20); it has %d children; branch points before it: %lld\n",
           (long long)target, (long long)index, children, (long long)branch_points_before);
    return 0;
}
