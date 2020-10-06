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
    /// (Y^2)Z=X^3+AXZ^2+BZ^3
    /// </summary>
    class PointProjective : IPoint
    {
        public BigInteger x { get; private set; }
        public BigInteger y { get; private set; }
        public BigInteger z { get; private set; }

        public EllipticCurveWeierstrass E { get; private set; }

        public PointProjective(BigInteger x, BigInteger y, BigInteger z, EllipticCurveWeierstrass E)
        {
            this.x = x;
            this.y = y;
            this.z = z;
            this.E = E;
        }

        public PointProjective(BigInteger x, BigInteger y, EllipticCurveWeierstrass E)
        {
            this.x = x;
            this.y = y;
            this.z = 1;
            this.E = E;
        }

        public static bool operator ==(PointProjective P, PointProjective Q)
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
        public static bool operator !=(PointProjective P, PointProjective Q)
        {
            return !(P == Q);
        }
        public static PointProjective operator -(PointProjective P)
        {
            return new PointProjective(P.x, -P.y, P.z, P.E as EllipticCurveWeierstrass);
        }
        public static PointProjective operator +(PointProjective P, PointProjective Q)//12M + 2S + 6add + 1*2
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
                BigInteger X3 =( r * r - J - 2 * V) % P.E.n;
                BigInteger Y3 = (r * (V - X3) - 2 * P.y * J) % P.E.n;
                BigInteger t = P.z + H;
                BigInteger Z3 =( 2 - Z1Z1 - HH) % P.E.n;
                return new PointProjective(X3, Y3, Z3, P.E as EllipticCurveWeierstrass);
            }
            else
            {
                BigInteger Y1Z2 = P.y * Q.z;
                BigInteger X1Z2 = P.x * Q.z;
                BigInteger Z1Z2 = P.z * Q.z;  //e
                BigInteger u = Q.y * P.z - Y1Z2;  //a
                BigInteger uu = u * u;
                BigInteger v = Q.x * P.z - X1Z2;  //b
                BigInteger vv = v * v;
                BigInteger vvv = v * vv;
                BigInteger R = vv * X1Z2;
                BigInteger A = uu * Z1Z2 - vvv - 2 * R;
                BigInteger x = (v * A) % P.E.n;
                BigInteger y = (u * (R - A) - vvv * Y1Z2) % P.E.n;
                BigInteger z = (vvv * Z1Z2) % P.E.n;
                return new PointProjective(x, y, z, P.E as EllipticCurveWeierstrass);
            }
        }
        public static PointProjective operator -(PointProjective P, PointProjective Q)
        {
            return P + (-Q);
        }
        public static PointProjective operator *(BigInteger k, PointProjective P)
        {
            PointProjective R = new PointProjective(0, 1, 0, P.E);
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

        public static PointProjective Double(PointProjective P)//5M + 6S + 1*a + 7add + 3*2 + 1*3
        {
            if (P.IsInfinite())
                return new PointProjective(0, 1, 0, P.E);
            BigInteger XX = P.x * P.x;
            BigInteger ZZ = P.z * P.z;
            BigInteger w = P.E.a * ZZ + 3 * XX;
            BigInteger s = 2 * P.y * P.z;
            BigInteger ss = s * s;
            BigInteger sss = s * ss;
            BigInteger R = P.y * s;
            BigInteger RR = R * R;
            BigInteger B = (P.x + R) * (P.x + R) - XX - RR;
            BigInteger h = w * w - 2 * B;
            BigInteger x = (h * s) % P.E.n;
            BigInteger y = (w * (B - h) - 2 * RR) % P.E.n;
            BigInteger z = sss % P.E.n;
            return new PointProjective(x, y, z, P.E);
        }
        public bool IsOnCurve()
        {
            return ((y * y * z - x * x * x - E.a * x * z * z - E.b * z * z * z) % E.n == 0);
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
            return n * (this as PointProjective);
        }
        public IPoint Add(IPoint P)
        {
            return this + (P as PointProjective);
        }
        public IPoint Double()
        {
            return this.Double();
        }
        public IPoint Clone()
        {
            return new PointProjective(this.x, this.y, this.z, this.E);
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
