using ECM2018.core.algorithms;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace ECM2018.core
{
    abstract class ECM
    {
        public static bool stop;
        protected static BigInteger StageOneAffine(ref IPoint P, BigInteger n, int B1)
        {
            var primes = Algorithms.Sieve(B1);
            for (int i = 0; i < primes.Length; i++)
            {
                BigInteger p = primes[i];
                BigInteger result = p;
                while (result < B1)
                {
                    P = P.Multiply(p);
                    if (P.z > 1)
                    {
                        return Algorithms.GCD(P.z, n);
                    }
                    if (P.IsInfinite())
                        return 0;
                    result *= p;
                    if (stop)
                        return 0;
                }
            }
            return 1;
        }
        protected static BigInteger StageOneProjective(ref IPoint P, BigInteger n, int B1)
        {
            var primes = Algorithms.Sieve(B1);
            int count = 0;
            BigInteger T = 1;
            for (int i = 0; i < primes.Length; i++)
            {
                BigInteger p = primes[i];
                BigInteger result = p;
                while (result < B1)
                {
                    P = P.Multiply(p);
                    T *= P.z;
                    T %= n;
                    if (count > 100)
                    {
                        if (T != 0)// вместо Т было P.z
                        {
                            BigInteger D = Algorithms.GCD(T, n);// вместо Т было P.z
                            if (D > 1)
                            {

                                return D;

                            }
                        }
                        else return 0;
                        count = 0;
                    }
                    if (P.IsInfinite())
                        return 0;
                    count++;
                    result *= p;
                    if (stop)
                        return 0;
                }
            }
            return 1; // результат не найден
        }
        protected static BigInteger StageOneEdwardsProjective(ref IPoint P, BigInteger n, int B1)
        {
            var primes = Algorithms.Sieve(B1);
            int count = 0;
            BigInteger T = 1;
            for (int i = 0; i < primes.Length; i++)
            {
                BigInteger p = primes[i];
                BigInteger result = p;
                while (result < B1)
                {
                    P = P.Multiply(p);
                    T *= P.x;
                    T %= n;
                    if (count > 100)
                    {
                        if (T != 0)// вместо Т было P.z
                        {
                            BigInteger D = Algorithms.GCD(T, n);// вместо Т было P.z
                            if (D > 1)
                                return D;
                        }
                        else return 0;
                        count = 0;
                    }
                    if (P.IsInfinite())
                        return 0;
                    count++;
                    result *= p;
                    if (stop)
                        return 0;
                }
            }
            return 1; // результат не найден
        }

        protected static BigInteger StageTwoAffine(ref IPoint P, BigInteger n, int B1, int B2)
        {
            var primes = Algorithms.Sieve(B1, B2);
            int maxDistance = 1;
            for (int i = 1; i < primes.Length; i++)
            {
                if (maxDistance < primes[i] - primes[i - 1])
                {
                    maxDistance = (int)(primes[i] - primes[i - 1]);
                }
            }
            var delta = new IPoint[maxDistance + 1];
            for (int i = 1; i < delta.Length; i++)
            {
                delta[i] = P.Multiply(i);
                if (delta[i].z > 1)
                {
                    return Algorithms.GCD(delta[i].z, n);
                }
            }

            IPoint S;
            S = P.Multiply(primes[0]);
            if (S.z > 1)
            {
                return Algorithms.GCD(S.z, n);
            }

            for (int i = 1; i < primes.Length; i++)
            {
                S = S.Add(delta[(int)(primes[i] - primes[i - 1])]);
                if (S.z > 1)
                {
                    return Algorithms.GCD(S.z, n);
                }
                if (S.IsInfinite())
                {
                    return 0;
                }
                if (stop)
                    return 0;
            }
            return 0;
        }
        protected static BigInteger StageTwoProjective(ref IPoint P, BigInteger n, int B1, int B2)
        {
            var primes = Algorithms.Sieve(B1, B2);
            int maxDistance = 1;
            for (int i = 1; i < primes.Length; i++)
            {
                if (maxDistance < primes[i] - primes[i - 1])
                {
                    maxDistance = (int)(primes[i] - primes[i - 1]);
                }
            }
            var delta = new IPoint[maxDistance + 1];
            delta[1] = P.Clone();
            delta[1].ToPointAffine();
            for (int i = 2; i < delta.Length; i++)
            {
                delta[i] = delta[i - 1].Add(delta[1]);
                //delta[i].ToPointAffine();

            }

            IPoint S;
            S = P.Multiply(primes[0]);
            BigInteger ProdZ = S.z;
            int count = 0;
            BigInteger d;
            for (int i = 1; i < primes.Length; i++)
            {
                S = S.Add(delta[(int)(primes[i] - primes[i - 1])]);
                ProdZ *= S.z;
                ProdZ %= n;
                if (count > 100)
                {
                    if (ProdZ == 0)
                    {
                        return 0;
                    }
                    d = Algorithms.GCD(ProdZ, n);
                    if (d > 1)
                    {
                        return d;
                    }
                    count = 0;
                }
                count++;
                if (stop)
                    return 0;
            }
            if (ProdZ == 0)
            {
                return 0;
            }
            else
            {
                d = Algorithms.GCD(ProdZ, n);
                if (d > 1)
                {
                    return d;
                }
            }
            return 0;
        }
        protected static BigInteger StageTwoProjective2(ref IPoint P, BigInteger n, int B1, int B2)
        {
            var primes = Algorithms.Sieve(B1, B2);
            BigInteger ProdZ = 1;
            int control = 0;
            for (int i = 0; i < primes.Length; i++)
            {
                P = P.Multiply(primes[i]);
                ProdZ *= P.z;
                ProdZ %= n;
                control++;
                if (control > 100)
                {
                    ProdZ = Algorithms.GCD(ProdZ, n);
                    if (ProdZ > 1 && ProdZ < n)
                    {
                        return ProdZ;
                    }
                    control = 0;
                }
            }
            ProdZ = Algorithms.GCD(ProdZ, n);
            if (ProdZ > 1 && ProdZ < n)
            {
                return ProdZ;
            }
            return 0;

        }
        protected static BigInteger StageTwoEdwardsProjective(ref IPoint P, BigInteger n, int B1, int B2)
        {
            var primes = Algorithms.Sieve(B1, B2);
            int maxDistance = 1;
            for (int i = 1; i < primes.Length; i++)
            {
                if (maxDistance < primes[i] - primes[i - 1])
                {
                    maxDistance = (int)(primes[i] - primes[i - 1]);
                }
            }
            var delta = new IPoint[maxDistance + 1];
            delta[1] = P.Clone();
            delta[1].ToPointAffine();
            for (int i = 2; i < delta.Length; i++)
            {
                delta[i] = delta[i - 1].Add(delta[1]);
                //delta[i].ToPointAffine();
            }


            IPoint S;
            S = P.Multiply(primes[0]);
            BigInteger ProdZ = S.x;
            int count = 0;
            BigInteger d;
            for (int i = 1; i < primes.Length; i++)
            {
                S = S.Add(delta[(int)(primes[i] - primes[i - 1])]);
                if (S.x!=0)
                {
                    ProdZ *= S.x;
                    ProdZ %= n;

                }
                if (count > 100)
                {

                    if (ProdZ == 0)
                    {
                        return 0;
                    }
                    d = Algorithms.GCD(ProdZ, n);
                    if (d > 1)
                    {
                            return d;
                    }
                    count = 0;
                }
                count++;
                if (stop)
                    return 0;
            }
            if (ProdZ == 0)
            {
                return 0;
            }
            else
            {
                d = Algorithms.GCD(ProdZ, n);
                if (d > 1)
                {
                    return d;
                }
            }
            return 0;
        }

    }
}
