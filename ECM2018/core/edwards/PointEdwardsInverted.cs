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
    /// (X^2+AY^2)Z^2=(X^2)(Y^2)+D(Z^4)
    /// </summary>
    class PointEdwardsInverted : IPoint
    {
        public BigInteger x { get; private set; }
        public BigInteger y { get; private set; }
        public BigInteger z { get; private set; }

        public EllipticCurveEdwards E { get; private set; }

        public PointEdwardsInverted(BigInteger x, BigInteger y, BigInteger z, EllipticCurveEdwards E)
        {
            this.x = x;
            this.y = y;
            this.z = z;
            this.E = E;
        }

        public PointEdwardsInverted(BigInteger x, BigInteger y, EllipticCurveEdwards E)
        {
            this.x = Algorithms.Invert(x, E.n);
            this.y = Algorithms.Invert(y, E.n);
            this.z = 1;
            this.E = E;
        }

        public static bool operator ==(PointEdwardsInverted P, PointEdwardsInverted Q)
        {
            return P.x == Q.x && P.y == Q.y && P.z == Q.z && P.E == Q.E;
        }
        public static bool operator !=(PointEdwardsInverted P, PointEdwardsInverted Q)
        {
            return !(P == Q);
        }
        public static PointEdwardsInverted operator -(PointEdwardsInverted P)
        {
            return new PointEdwardsInverted(-P.x, P.y, P.z, P.E);
        }
        public static PointEdwardsInverted operator +(PointEdwardsInverted P, PointEdwardsInverted Q)//9M + 1S + 1*a + 1*d + 7add
        {

            BigInteger A, B, C, D, E, F, G, H, I, x, y, z;
            if (Q.z == 1)
            {
                B = P.E.b * P.z * P.z;
                C = P.x * Q.x;
                D = P.y * Q.y;
                E = C * D;
                H = C - P.E.a * D;
                I = (P.x + P.y) * (Q.x + Q.y) - C - D;
                x = ((E + B) * H) % P.E.n;
                y = ((E - B) * I) % P.E.n;
                z = (P.z * H * I) % P.E.n;
                return new PointEdwardsInverted(x, y, z, P.E);
            }
            A = P.z * Q.z;
            B = P.E.b * A * A;
            C = P.x * Q.x;
            D = P.y * Q.y;
            E = C * D;
            H = C - P.E.a * D;
            I = (P.x + P.y) * (Q.x + Q.y) - C - D;
            x = ((E + B) * H) % P.E.n;
            y = ((E - B) * I) % P.E.n;
            z = (A * H * I) % P.E.n;
            return new PointEdwardsInverted(x, y, z, P.E);
        }
        public static PointEdwardsInverted operator -(PointEdwardsInverted P, PointEdwardsInverted Q)
        {
            return P + (-Q);
        }
        public static PointEdwardsInverted operator *(BigInteger k, PointEdwardsInverted P)
        {
            bool isFirst = true;
            PointEdwardsInverted R = new PointEdwardsInverted(0, 1, 1, P.E);

            while (k > 0)
            {
                if (k % 2 == 1 && isFirst)
                {
                    R = new PointEdwardsInverted(P.x, P.y, P.z, P.E);
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

        public static PointEdwardsInverted Double(PointEdwardsInverted P)//3M + 4S + 1*a + 1*d2 + 6add
        {

            BigInteger A = P.x * P.x;
            BigInteger B = P.y * P.y;
            BigInteger U = P.E.a * B;
            BigInteger C = A + U;
            BigInteger D = A - U;
            BigInteger W = P.x + P.y;
            BigInteger E = W * W - A - B;
            BigInteger x = (C * D) % P.E.n;
            BigInteger y = (E * (C - 2 * P.E.b * P.z * P.z)) % P.E.n;
            BigInteger z = (D * E) % P.E.n;
            return new PointEdwardsInverted(x, y, z, P.E);
        }
        public bool IsOnCurve()
        {
            return (((E.a * y * y + x * x) * z * z - x * x * y * y - E.b * z * z * z * z) % E.n == 0);
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
            return n * (this as PointEdwardsInverted);
        }
        public IPoint Add(IPoint P)
        {
            return this + (P as PointEdwardsInverted);
        }
        public IPoint Double()
        {
            return this.Double();
        }
        public IPoint Clone()
        {
            return new PointEdwardsInverted(this.x, this.y, this.z, this.E);
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
