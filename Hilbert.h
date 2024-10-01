#pragma once

#if 0

#include <vector>

// TODO: look for todos below
// TODO: does anything need cleaning up?
// TODO: should we make these functions static cause they don't affect storage?
// TODO: maybe make a new repo for this 2d code and tests

// Adapted from "Compact indices" section of https://pdebuyl.be/blog/2015/hilbert-curve.html

template <int N>
class Hilbert
{
public:
    // TODO: where should this go?
    void Test()
    {
        int indexCount = 0;
        for (int i : compact_M)
            indexCount += i;
        indexCount = 1 << indexCount;

        for (int index = 0; index < indexCount; ++index)
        {
            int point[3];
            IndexToPoint(index, point);
            int indexRT = PointToIndex(point);
            int ijkl = 0;
        }
    }
    
    // Compute the compact Hilbert index for point p
    // TR_algo7 in the python.
    int PointToIndex(const int P[N])
    {
        int h = 0;
        int ve = 0;
        int vd = 2;
        int m = 0;
        for (int i : compact_M)
            m = std::max(m, i);
        for (int i = m - 1; i >= 0; i--)
        {
            int mu = extract_mask(i);
            
            int mu_norm = 0;
            for (int j = 0; j < N; ++j)
                mu_norm += bit_component(mu, j);

            mu = rotate_right(mu, vd + 1);

            int pi = rotate_right(ve, vd + 1) & ((~mu) & (1 << N) - 1);

            // 2. construct a integer whose bits are given by l
            int l = 0;
            for (int j = 0; j < N; ++j)
            {
                if (bit_component(P[j], i))
                    l |= (1 << j);
            }

            l = T(ve, vd, l);
            int w = inverse_gc(l);

            int r = gcr(w, mu, pi);

            ve = ve ^ rotate_left(e(w), vd + 1);

            vd = (vd + d(w) + 1) % N;

            h = (h << mu_norm) | r;
        }
        return h;
    }

    // Compute the point with compact Hilbert index h
    // TR_algo8 in the python.
    void IndexToPoint(int index, int P[N])
    {
        int h = 0; // TODO: i had to add h!
        int ve = 0;
        int vd = 2;
        int k = 0;
        for (int i = 0; i < N; ++i)
            P[i] = 0;

        int m = 0;
        int vM = 0;
        for (int i : compact_M)
        {
            m = std::max(m, i);
            vM += i;
        }

        for (int i = m - 1; i >= 0; i--)
        {
            int mu = extract_mask(i);

            int mu_norm = 0;
            for (int j = 0; j < N; ++j)
                mu_norm += bit_component(mu, j);

            mu = rotate_right(mu, vd + 1);

            int pi = rotate_right(ve, vd + 1) & (~mu & (1 << N) - 1);

            int r = 0;
            for (int j = 0; j < mu_norm; ++j)
            {
                if (bit_component(h, vM - k - (j + 1))) // TODO: where does h come from?
                    r |= (1 << j);
            }
            k = k + mu_norm;
            int w = gcr_inv(r, mu, pi);
            int l = gc(w);
            l = T_inv(ve, vd, l);
            for (int j = 0; j < N; ++j)
                P[j] |= (bit_component(l, j) << i);

            ve = ve ^ (rotate_left(e(w), vd + 1));
            vd = (vd + d(w) + 1) % N;
        }
    }

private:
    // TODO: i think this needs to be a parameter. maybe variadic template arguments, where the length defines N
    // # definition of the size of the space. compact_M is of length N
    static inline const int compact_M[3] = { 3, 2, 2 };

    // Return the Gray code index of i.
    int gc(int i)
    {
        return i ^ (i >> 1);
    }

    // The inverse gray code.
    int inverse_gc(int g)
    {
        int i = g;
        int j = 1;
        while (j < N)
        {
            i = i ^ (g >> j);
            j = j + 1;
        }
        return i;
    }

    // The direction between subcube i and the next one
    int g(int i)
    {
        return int(std::log2(gc(i) ^ gc(i + 1)));
    }

    // The direction of the arrow whithin a subcube.
    int d(int i)
    {
        if (i == 0)
            return 0;
        else if ((i % 2) == 0)
            return g(i - 1) % N;
        else
            return g(i) % N;
    }

    //Transform b.
    int T(int e, int d, int b)
    {
        int out = b ^ e;
        return rotate_right(out, d + 1);
    }

    // Inverse transform b.
    int T_inv(int e, int d, int b)
    {
        return T(rotate_right(e, d + 1), N - d - 2, b);
    }

    // Return the entry point of hypercube i.
    int e(int i)
    {
        if (i == 0)
            return 0;
        else
            return gc(2 * int((i - 1) / 2)); // TODO: i think this is just gc(i-1)? after things are working try that cleanup?
    }

    // Compute the gray code rank of i given the mask mu.
    // Algorithm 4 in [TR]
    int gcr(int i, int mu, int pi)
    {
        int r = 0;
        for (int k = N - 1; k >= 0; k--)
        {
            if (bit_component(mu, k))
                r = (r << 1) | bit_component(i, k);
        }
        return r;
    }

    // Inverse of the gray code rank, given the mask mu and the pattern pi.
    // Algorithm 5 in [TR]
    int gcr_inv(int r, int mu, int pi)
    {
        int i = 0;
        int g = 0;
        int j = -1;
        for (int k = 0; k < N; ++k)
            j += bit_component(mu, k);
        for (int k = N - 1; k >= 0; k--)
        {
            if (bit_component(mu, k))
            {
                i |= (bit_component(r, j) << k);
                g |= (((bit_component(i, k) + bit_component(i, k + 1)) % 2) << k);
                j -= 1;
            }
            else
            {
                g |= (bit_component(pi, k) << k);
                i |= (((bit_component(g, k) + bit_component(i, k + 1)) % 2) << k);
            }
        }
        return i;
    }

    // Extract the mask for iteration i of the algorithm.
    // Algorithm 6 in [TR]
    int extract_mask(int i)
    {
        int mu = 0;
        for (int j = N - 1; j >= 0; j--)
        {
            mu = mu << 1;
            if (compact_M[j] > i)
                mu = mu | 1;
        }
        return mu;
    }

    // Rotate x by d bits to the right.
    int rotate_right(int x, int d)
    {
        d = d % N;
        int out = x >> d;
        for (int i = 0; i < d; i++)
        {
            int bit = (x & (1 << i)) >> i;
            out |= bit << (N + i - d);
        }
        return out;
    }

    // Rotate x by d bits to the left.
    int rotate_left(int x, int d)
    {
        d = d % N;
        int out = x << d;
        int excess = out;
        out = out & ((1 << N) - 1);
        for (int i = 0; i < d; ++i)
        {
            int bit = (x & (1 << (N - 1 - d + 1 + i))) >> (N - 1 - d + 1 + i);
            out |= bit << i;
        }
        return out;
    }

    // Return i-th bit of x
    int bit_component(int x, int i)
    {
        return (x & (1 << i)) ? 1 : 0;
    }
};

#endif

// From https://hugocisneros.com/notes/hilbert_curve_indexing/

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
