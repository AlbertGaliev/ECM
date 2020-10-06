using ECM2018.core.algorithms;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace ECM2018.core.edwards
{
    /// <summary>
    /// AX^2+Y^2=Z^2+D(T^2)
    /// </summary>
    class PointEdwardsExtended : IPoint
    {
        public BigInteger x { get; private set; }
        public BigInteger y { get; private set; }
        public BigInteger z { get; private set; }
        public BigInteger t { get; private set; }

        public EllipticCurveEdwards E { get; private set; }

        public PointEdwardsExtended(BigInteger x, BigInteger y, BigInteger z, BigInteger t, EllipticCurveEdwards E)
        {
            this.x = x;
            this.y = y;
            this.z = z;
            this.t = t;
            this.E = E;
        }

        public PointEdwardsExtended(BigInteger x, BigInteger y, EllipticCurveEdwards E)
        {
            this.x = x;
            this.y = y;
            this.z = 1;
            this.t = (x * y) % E.n;
            this.E = E;
        }

        public static bool operator ==(PointEdwardsExtended P, PointEdwardsExtended Q)
        {
            return P.x == Q.x && P.y == Q.y && P.z == Q.z && P.E == Q.E;
        }
        public static bool operator !=(PointEdwardsExtended P, PointEdwardsExtended Q)
        {
            return !(P == Q);
        }
        public static PointEdwardsExtended operator -(PointEdwardsExtended P)
        {
            return new PointEdwardsExtended(-P.x, P.y, P.z, -P.t, P.E);
        }
        public static PointEdwardsExtended operator +(PointEdwardsExtended P, PointEdwardsExtended Q)//9M + 1*a + 7add
        {

            BigInteger A, B, C, D, E, F, G, H, I, x, y, z, t;
            if (Q.z == 1)
            {
                A = P.x * Q.x;
                B = P.y * Q.y;
                C = P.z * Q.t;
                D = P.t;
                E = D + C;
                F = (P.x - P.y) * (Q.x + Q.y) + B - A;
                G = B + P.E.a * A;
                H = D - C;
                x = (E * F) % P.E.n;
                y = (G * H) % P.E.n;
                t = (E * H) % P.E.n;
                z = (F * G) % P.E.n;
                return new PointEdwardsExtended(x, y, z, t, P.E);
            }
            A = P.x * Q.x;
             B = P.y * Q.y;
             C = P.z * Q.t;
             D = P.t * Q.z;
             E = D + C;
             F = (P.x - P.y) * (Q.x + Q.y) + B - A;
             G = B + P.E.a * A;
             H = D - C;
             x = (E * F) % P.E.n;
             y = (G * H) % P.E.n;
             t = (E * H) % P.E.n;
             z = (F * G) % P.E.n;
            return new PointEdwardsExtended(x, y, z, t, P.E);
        }
        public static PointEdwardsExtended operator -(PointEdwardsExtended P, PointEdwardsExtended Q)
        {
            return P + (-Q);
        }
        public static PointEdwardsExtended operator *(BigInteger k, PointEdwardsExtended P)
        {
            var R = new PointEdwardsExtended(0, 1, P.E);
            if (P.IsInfinite())
                return R;
            while (k > 0)
            {
                if (k % 2 == 1)
                    R += P;
                k /= 2;
                P = Double(P);
            }
            return R;
        }

        public static PointEdwardsExtended Double(PointEdwardsExtended P)//4M + 4S + 1*a + 6add + 1*2
        {
            BigInteger A = P.x * P.x;
            BigInteger B = P.y * P.y;
            BigInteger C = 2 * P.z * P.z;
            BigInteger D = P.E.a * A;
            BigInteger W = P.x + P.y;
            BigInteger E = W * W - A - B;
            BigInteger G = D + B;
            BigInteger F = G - C;
            BigInteger H = D - B;
            BigInteger x = (E * F) % P.E.n;
            BigInteger y = (G * H) % P.E.n;
            BigInteger t = (E * H) % P.E.n;
            BigInteger z = (F * G) % P.E.n;
            return new PointEdwardsExtended(x, y, z, t, P.E);
        }

        public static PointEdwardsExtended Tripple(PointEdwardsExtended P)//4M + 4S + 1*a + 6add + 1*2
        {
            BigInteger YY = P.y * P.y;
            BigInteger aXX = P.E.a * P.x * P.y;
            BigInteger Ap = YY + aXX;
            BigInteger B = 2 * (2 * P.z * P.z - Ap);
            BigInteger xB = aXX * B;
            BigInteger yB = YY * B;
            BigInteger AA = Ap * (YY - aXX);
            BigInteger F = AA - yB;
            BigInteger G = AA + xB;
            BigInteger xE = P.x * (yB + AA);
            BigInteger yH = P.y * (xB - AA);
            BigInteger zF = P.z * F;
            BigInteger zG = P.z * G;
            BigInteger x = (xE * zF) % P.E.n;
            BigInteger y = (yH * zG) % P.E.n;
            BigInteger z = (zF * zG) % P.E.n;
            BigInteger t = (xE * yH) % P.E.n;
            return new PointEdwardsExtended(x, y, z, t, P.E);
        }
        public bool IsOnCurve()
        {
            return ((E.a * x * x + y * y - z * z - E.b * t * t) % E.n == 0);
        }
        public bool IsInfinite()
        {
            return (x == 0 && y != 0 && z == 0) || (x != 0 && y == 0 && z == 0);
        }
        public override bool Equals(Object obj)
        {
            throw new NotImplementedException();
        }
        public override int GetHashCode()
        {
            return this.GetHashCode();
        }
        public IPoint Multiply(BigInteger n)
        {
            return n * (this as PointEdwardsExtended);
        }
        public IPoint Add(IPoint P)
        {
            return this + (P as PointEdwardsExtended);
        }
        public IPoint Double()
        {
            return this.Double();
        }
        public IPoint Clone()
        {
            return new PointEdwardsExtended(this.x, this.y, this.z, this.t, this.E);
        }
        public void ToPointAffine()
        {
            if (this.z != 0)
            {
                BigInteger inv = Algorithms.Invert(this.z, this.E.n);
                this.x = (this.x * inv) % this.E.n;
                this.y = (this.y * inv) % this.E.n;
                this.z = 1;
                this.z = (x*y) % this.E.n;
            }
        }

    }

}
