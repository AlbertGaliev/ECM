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
    /// y^2=x^3+ax+b
    /// </summary>
    class PointAffine : IPoint
    {
        public BigInteger x { get; private set; }
        public BigInteger y { get; private set; }
        public BigInteger z { get; private set; }

        public EllipticCurveWeierstrass E { get; private set; }

        public PointAffine(BigInteger x, BigInteger y, BigInteger z, EllipticCurveWeierstrass E)
        {
            this.x = x;
            this.y = y;
            this.z = z;
            this.E = E;
        }

        public PointAffine(BigInteger x, BigInteger y, EllipticCurveWeierstrass E)
        {
            this.x = x;
            this.y = y;
            this.z = 1;
            this.E = E;
        }

        public static bool operator ==(PointAffine P, PointAffine Q)
        {
            return P.x == Q.x && P.y == Q.y && P.z == Q.z && P.E == Q.E;
        }
        public static bool operator !=(PointAffine P, PointAffine Q)
        {
            return !(P == Q);
        }
        public static PointAffine operator -(PointAffine P)
        {
            return new PointAffine(P.x, -P.y, P.z, P.E);
        }
        public static PointAffine operator +(PointAffine P, PointAffine Q)
        {
            BigInteger temp;
            BigInteger lambda;
            BigInteger inv;
            BigInteger d;
            BigInteger x;
            BigInteger y;
            
            if (P.x < 0)
                P.x += P.E.n;
            if (Q.x < 0)
                Q.x += Q.E.n;

            if (P.IsInfinite())
                return Q;
            else if (Q.IsInfinite())
                return P;

            else if (P == -Q || (P == Q && P.y == 0))
                return new PointAffine(0, 1, 0, P.E);
            else if (P == Q)
                return Double(P);
            else if (P.x == Q.x)    // случай, если у-координата противоположной точки неотрицательна
            {
                temp = (P.y + Q.y) % P.E.n;
                if (temp == 0)
                    return new PointAffine(0, 1, 0, P.E);
                if (temp < 0)  
                    temp += P.E.n;
                Algorithms.ExtendedEuclidAlgorithm(temp, P.E.n, out inv, out d);
                if (d > 1 && d < P.E.n)
                    return new PointAffine(0, 0, temp, P.E);    //делитель найден, можно прекращать вычисления
                lambda = (3 * P.x * P.x + P.E.a) * inv;
            }
            else
            {
                temp = (Q.x - P.x) % P.E.n;
                if (temp < 0)
                    temp += P.E.n;
                Algorithms.ExtendedEuclidAlgorithm(temp, P.E.n, out inv, out d);
                if (d > 1 && d < P.E.n)
                    return new PointAffine(0, 0, temp, P.E);    //делитель найден, можно прекращать вычисления
                lambda = (Q.y - P.y) * inv;
            }
            x = (lambda * lambda - P.x - Q.x) % P.E.n;
            y = (lambda * (P.x - x) - P.y) % P.E.n;
            return new PointAffine(x, y, 1, P.E);
        }
        public static PointAffine operator -(PointAffine P, PointAffine Q)
        {
            return P + (-Q);
        }
        public static PointAffine operator *(BigInteger k, PointAffine P)
        {

            PointAffine R = new PointAffine(0, 1, 0, P.E);
            if (P.IsInfinite())
                return R;
            while (k > 0)
            {
                if (P.z > 1)    //делитель найден, можно прекращать вычисления
                    return P;
                if (k % 2 == 1)
                    R = P + R;
                k /= 2;
                P = Double(P);//P = P + P;
            }
            return R;
        }

        public static PointAffine Double(PointAffine P)
        {
            if (P.y == 0)
                return new PointAffine(0, 1, 0, P.E);
            BigInteger temp = (2 * P.y) % P.E.n;
            if (temp < 0)
                temp += P.E.n;
            BigInteger inv;
            BigInteger d;
            Algorithms.ExtendedEuclidAlgorithm(temp, P.E.n, out inv, out d);
            if (d > 1 && d < P.E.n)
                return new PointAffine(0, 0, temp, P.E);    //делитель найден, можно прекращать вычисления
            BigInteger lambda = (3 * P.x * P.x + P.E.a) * inv;
            BigInteger x = (lambda * lambda - 2 * P.x) % P.E.n;
            BigInteger y = (lambda * (P.x - x) - P.y) % P.E.n;
            return new PointAffine(x, y, 1, P.E);
        }
        public bool IsOnCurve()
        {
            if (IsInfinite())
                return true;
            return ((y * y - x * x * x - E.a * x - E.b) % E.n == 0);
        }
        public bool IsInfinite()
        {
            return (x == 0 && y == 1 && z == 0);
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
            return n * (this as PointAffine);
        }
        public IPoint Add(IPoint P)
        {
            return this + (P as PointAffine);
        }
        public IPoint Double()
        {
            return this.Double();
        }
        public IPoint Clone()
        {
            return new PointAffine(this.x, this.y, this.z, this.E);
        }

        
        public void ToPointAffine()
        {

        }

    }

}
