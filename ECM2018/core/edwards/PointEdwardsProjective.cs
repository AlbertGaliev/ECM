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
    /// (AX^2+Y^2)Z^2=Z^4+D(X^2)(Y^2)
    /// </summary>
    class PointEdwardsProjective : IPoint
    {
        public BigInteger x { get; private set; }
        public BigInteger y { get; private set; }
        public BigInteger z { get; private set; }

        public EllipticCurveEdwards E { get; private set; }

        public PointEdwardsProjective(BigInteger x, BigInteger y, BigInteger z, EllipticCurveEdwards E)
        {
            this.x = x;
            this.y = y;
            this.z = z;
            this.E = E;
        }

        public PointEdwardsProjective(BigInteger x, BigInteger y, EllipticCurveEdwards E)
        {
            this.x = x;
            this.y = y;
            this.z = 1;
            this.E = E;
        }

        public static bool operator ==(PointEdwardsProjective P, PointEdwardsProjective Q)
        {
            if (P.E == Q.E)
            {
                BigInteger u1 = (P.y * Q.z) % P.E.n;
                BigInteger u2 = (P.z * Q.y) % P.E.n;

                BigInteger v1 = (P.x * Q.z) % P.E.n;
                BigInteger v2 = (P.z * Q.x) % P.E.n;

                if (u1 < 0) u1 += P.E.n;
                if (u2 < 0) u2 += P.E.n;
                if (v1 < 0) v1 += P.E.n;
                if (v2 < 0) v2 += P.E.n;

                if (u1 == u2)
                    if (v1 == v2)
                        return true;
            }
            return false;
        }
        public static bool operator !=(PointEdwardsProjective P, PointEdwardsProjective Q)
        {
            return !(P == Q);
        }
        public static PointEdwardsProjective operator -(PointEdwardsProjective P)
        {
            return new PointEdwardsProjective(-P.x, P.y, P.z, P.E);
        }
        public static PointEdwardsProjective operator +(PointEdwardsProjective P, PointEdwardsProjective Q)//10M + 1S + 1*a + 1*d + 7add
        {
            BigInteger A, B, C, D, E, F, G, x, y, z;
            if (Q.z == 1)
            {
                B = P.z * P.z;
                C = P.x * Q.x;
                D = P.y * Q.y;
                E = P.E.b * C * D;
                F = B - E;
                G = B + E;
                x = (P.z * F * ((P.x + P.y) * (Q.x + Q.y) - C - D)) % P.E.n;
                y = (P.z * G * (D - P.E.a * C)) % P.E.n;
                z = (F * G) % P.E.n;
                return new PointEdwardsProjective(x, y, z, P.E);
            }
             A = P.z * Q.z;
             B = A * A;
             C = P.x * Q.x;
             D = P.y * Q.y;
             E = P.E.b * C * D;
             F = B - E;
             G = B + E;
             x = (A * F * ((P.x + P.y) * (Q.x + Q.y) - C - D)) % P.E.n;
             y = (A * G * (D - P.E.a * C)) % P.E.n;
             z = (F * G) % P.E.n;
            return new PointEdwardsProjective(x, y, z, P.E);
        }
        public static PointEdwardsProjective operator -(PointEdwardsProjective P, PointEdwardsProjective Q)
        {
            return P + (-Q);
        }
        public static PointEdwardsProjective operator *(BigInteger k, PointEdwardsProjective P)
        {
            bool isFirst = true;
            var R = new PointEdwardsProjective(0, 1, 1, P.E);

            while (k > 0)
            {
                if (k % 2 == 1 && isFirst)
                {
                    R = new PointEdwardsProjective(P.x, P.y, P.z, P.E);
                    isFirst = false;
                }
                else if (k % 2 == 1)
                {
                    R += P;
                }
                k /= 2;
                P = Double(P);
            }
            return R;
        }

        public static PointEdwardsProjective Double(PointEdwardsProjective P)//3M + 4S + 1*a + 6add + 1*2
        {
            BigInteger n = P.E.n;
            BigInteger A = (P.x + P.y);
            BigInteger B = A * A;
            BigInteger C = P.x * P.x;
            BigInteger D = P.y * P.y;
            BigInteger E = P.E.a * C;
            BigInteger F = E + D;
            BigInteger H = P.z * P.z;
            BigInteger J = F - 2 * H;
            BigInteger x = ((B - C - D) * J) % P.E.n;
            BigInteger y = (F * (E - D)) % P.E.n;
            BigInteger z = (F * J) % P.E.n;
            return new PointEdwardsProjective(x, y, z, P.E);
        }

        public static PointEdwardsProjective Tripple(PointEdwardsProjective P)//9M + 3S + 1*a + 7add + 2*2
        {
            BigInteger YY = P.y * P.y;
            BigInteger aXX = P.E.a * P.x * P.x;
            BigInteger Ap = YY + aXX;
            BigInteger B = 2 * (2 * P.z * P.z - Ap);
            BigInteger xB = aXX * B;
            BigInteger yB = YY * B;
            BigInteger AA = Ap * (YY - aXX);
            BigInteger F = AA - yB;
            BigInteger G = AA + xB;
            BigInteger x = (P.x * (yB + AA) * F) % P.E.n;
            BigInteger y = (P.y * (xB - AA) * G) % P.E.n;
            BigInteger z = (P.z * F * G) % P.E.n;
            return new PointEdwardsProjective(x, y, z, P.E);
        }
        public bool IsOnCurve()
        {
            return (((y * y + E.a * x * x) * z * z - z * z * z * z - E.b * x * x * y * y) % E.n == 0);
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
            return n * (this as PointEdwardsProjective);
        }
        public IPoint Add(IPoint P)
        {
            return this + (P as PointEdwardsProjective);
        }
        public IPoint Double()
        {
            return this.Double();
        }
        public IPoint Clone()
        {
            return new PointEdwardsProjective(this.x, this.y, this.z, this.E);
        }
        public void ToPointAffine()
        {
            if (this.z != 0)
            {
                BigInteger inv = Algorithms.Invert(this.z, this.E.n);
                this.x = (this.x * inv) % this.E.n;
                this.y = (this.y * inv) % this.E.n;
                this.z = 1;
            }
        }

    }

}

