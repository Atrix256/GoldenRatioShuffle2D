#pragma once

// From https://hugocisneros.com/notes/hilbert_curve_indexing/

namespace Hilbert
{
    //rotate/flip a quadrant appropriately
    static void rot(int n, int* x, int* y, int rx, int ry) {
        if (ry == 0) {
            if (rx == 1) {
                *x = n - 1 - *x;
                *y = n - 1 - *y;
            }

            //Swap x and y
            int t = *x;
            *x = *y;
            *y = t;
        }
    }

    //convert (x,y) to d
    static int xy2d(int n, int x, int y) {
        int rx, ry, s, d = 0;
        for (s = n / 2; s > 0; s /= 2) {
            rx = (x & s) > 0;
            ry = (y & s) > 0;
            d += s * s * ((3 * rx) ^ ry);
            rot(n, &x, &y, rx, ry);
        }
        return d;
    }

    //convert d to (x,y)
    static void d2xy(int n, int d, int* x, int* y) {
        int rx, ry, s, t = d;
        *x = *y = 0;
        for (s = 1; s < n; s *= 2) {
            rx = 1 & (t / 2);
            ry = 1 & (t ^ rx);
            rot(s, x, y, rx, ry);
            *x += s * rx;
            *y += s * ry;
            t /= 4;
        }
    }

};
