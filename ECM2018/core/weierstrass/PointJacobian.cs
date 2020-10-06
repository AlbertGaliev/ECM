using ECM2018.core.algorithms;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace ECM2018.core.weierstrass
{
    /// <summary>
    /// Y^2=X^3+AXZ^4+BZ^6
    /// </summary>
    class PointJacobian : IPoint
    {
        public BigInteger x { get; private set; }
        public BigInteger y { get; private set; }
        public BigInteger z { get; private set; }

        public EllipticCurveWeierstrass E { get; private set; }

        public PointJacobian(BigInteger x, BigInteger y, BigInteger z, EllipticCurveWeierstrass E)
        {
            this.x = x;
            this.y = y;
            this.z = z;
            this.E = E;
        }

        public PointJacobian(BigInteger x, BigInteger y, EllipticCurveWeierstrass E)
        {
            this.x = x;
            this.y = y;
            this.z = 1;
            this.E = E;
        }

        public static bool operator ==(PointJacobian P, PointJacobian Q)
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
        public static bool operator !=(PointJacobian P, PointJacobian Q)
        {
            return !(P == Q);
        }
        public static PointJacobian operator -(PointJacobian P)
        {
            return new PointJacobian(P.x, -P.y, P.z, P.E as EllipticCurveWeierstrass);
        }
        public static PointJacobian operator +(PointJacobian P, PointJacobian Q)//12M + 4S + 6add + 1*2
        {
            if (P.IsInfinite())
                return Q;
            if (Q.IsInfinite())
                return P;
            if (Q.z == 1)
            {
                BigInteger Z1Z1 = P.z * P.z;
                BigInteger U2 = Q.x * Z1Z1;
                BigInteger S2 = Q.y * P.z * Z1Z1;
                BigInteger H = U2 - P.x;
                BigInteger HH = H * H;
                BigInteger I = 4 * HH;
                BigInteger J = H * I;
                BigInteger r = 2 * (S2 - P.y);
                BigInteger V = P.x * I;
                BigInteger x = (r * r - J - 2 * V) % P.E.n;
                BigInteger y = (r * (V - x) - 2 * P.y * J) % P.E.n;
                BigInteger t = (P.z + H);
                BigInteger z = (t * t - Z1Z1 - HH) % P.E.n;
                return new PointJacobian(x, y, z, P.E as EllipticCurveWeierstrass);
            }
            else
            {
                BigInteger Z1Z1 = P.z * P.z;
                BigInteger Z2Z2 = Q.z * Q.z;
                BigInteger U1 = P.x * Z2Z2;
                BigInteger U2 = Q.x * Z1Z1;
                BigInteger S1 = P.y * Q.z * Z2Z2;
                BigInteger S2 = Q.y * P.z * Z1Z1;
                BigInteger H = U2 - U1;
                BigInteger HH = H * H;
                BigInteger HHH = H * HH;
                BigInteger r = S2 - S1;
                BigInteger V = U1 * HH;
                BigInteger x = (r * r - HHH - 2 * V) % P.E.n;
                BigInteger y = (r * (V - x) - S1 * HHH) % P.E.n;
                BigInteger z = (P.z * Q.z * H) % P.E.n;
                return new PointJacobian(x, y, z, P.E as EllipticCurveWeierstrass);
            }
        }
        public static PointJacobian operator -(PointJacobian P, PointJacobian Q)
        {
            return P + (-Q);
        }
        public static PointJacobian operator *(BigInteger k, PointJacobian P)
        {
            PointJacobian R = new PointJacobian(0, 1, 0, P.E);
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

        public static PointJacobian Double(PointJacobian P)//3M + 6S + 1*a + 4add + 2*2 + 1*3 + 1*4 + 1*8
        {
            if (P.IsInfinite())
                return new PointJacobian(0, 1, 0, P.E);
            BigInteger XX = P.x * P.x;
            BigInteger YY = P.y * P.y;
            BigInteger ZZ = P.z * P.z;
            BigInteger S = 4 * P.x * YY;
            BigInteger M = 3 * XX + P.E.a * ZZ * ZZ;
            BigInteger T = M * M - 2 * S;
            BigInteger x = T % P.E.n;
            BigInteger y = (M * (S - T) - 8 * YY * YY) % P.E.n;
            BigInteger z = (2 * P.y * P.z) % P.E.n;
            return new PointJacobian(x, y, z, P.E);
        }
        public bool IsOnCurve()
        {
            return ((y * y - x * x * x - E.a * x * z * z * z * z - E.b * z * z * z * z * z * z) % E.n == 0);
        }
        public bool IsInfinite()
        {
            return (x == 0 && z == 0 && y != 0);
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
            return n * (this as PointJacobian);
        }
        public IPoint Add(IPoint P)
        {
            return this + (P as PointJacobian);
        }
        public IPoint Double()
        {
            return this.Double();
        }
        public IPoint Clone()
        {
            return new PointJacobian(this.x, this.y, this.z, this.E);
        }

        public void ToPointAffine()
        {
            if (this.z != 0)
            {
                BigInteger inv = Algorithms.Invert(this.z, this.E.n);
                this.x = (this.x * inv * inv) % this.E.n;
                this.y = (this.y * inv * inv * inv) % this.E.n;
                this.z = 1;
            }
        }

    }
}
