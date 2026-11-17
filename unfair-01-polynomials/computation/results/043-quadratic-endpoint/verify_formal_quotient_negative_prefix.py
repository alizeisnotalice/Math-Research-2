#!/usr/bin/env python3
"""Integer formal-series regression for the positive binary calibration example."""

def quotient_series(numer, denom, terms):
    # denom[0]=1.  Solve denom * out = numer coefficient by coefficient.
    out = []
    for k in range(terms):
        target = numer[k] if k < len(numer) else 0
        value = target - sum(denom[j] * out[k - j]
                             for j in range(1, min(k, len(denom) - 1) + 1))
        out.append(value)
    return out


def main():
    # G/Q=(1+x)/(1+x^2), Q/G=(1+x^2)/(1+x).
    forward = quotient_series([1, 1], [1, 0, 1], 7)
    reverse = quotient_series([1, 0, 1], [1, 1], 7)
    assert forward == [1, 1, -1, -1, 1, 1, -1], forward
    assert reverse == [1, -1, 2, -2, 2, -2, 2], reverse
    # n=2, b_1=0; m=1, a_0=a_1=1.
    assert 1 + forward[2] + forward[1] * 0 == 0
    assert 1 + reverse[1] == 0
    print("G/Q coefficients through x^6:", forward)
    print("Q/G coefficients through x^6:", reverse)
    print("forward n=2 row and reverse m=1 row both equal zero")


if __name__ == "__main__":
    main()
