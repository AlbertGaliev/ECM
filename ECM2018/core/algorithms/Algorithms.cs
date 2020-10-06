using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace ECM2018.core.algorithms
{
    class Algorithms
    {
        static Random rand = new Random();
        public static BigInteger GCD(BigInteger a, BigInteger b)
        {
            if (a<0)
            {
                a += b;
            }
            if (a < b)
            {
                BigInteger temp = b;
                b = a;
                a = temp;
            }
            while (a % b != 0)
            {
                BigInteger temp = a % b;
                a = b;
                b = temp;
            }
            return BigInteger.Abs(b);
        }
        public static void ExtendedEuclidAlgorithm(BigInteger a, BigInteger b, out BigInteger inv, out BigInteger d)
        {
            BigInteger q, r, x1, x2, y1, y2, x, y;
            if (b == 0)
            {
                x = 1;
                y = 0;
                d = a;
                inv = x;
                return;
            }
            x1 = 0;
            x2 = 1;
            y1 = 1;
            y2 = 0;
            while (b > 0)
            {
                q = a / b;
                r = a - q * b;
                x = x2 - q * x1;
                y = y2 - q * y1;
                a = b;
                b = r;
                x2 = x1;
                x1 = x;
                y2 = y1;
                y1 = y;
            }
            x = x2;
            y = y2;
            d = a;
            inv = x;
        }

        public static bool MillerRabinTest(BigInteger R)
        {
            BigInteger k = (BigInteger)BigInteger.Log(R, 2);
            if (R == 2)
            {
                return true;
            }
            BigInteger t = R - 1;
            int s = 0;
            while (t % 2 == 0)
            {
                s++;
                t = t / 2;
            }
            for (BigInteger i = 0; i < k; i++)
            {
                BigInteger a = RandomBigInteger(R - 1);
                if (R % a == 0)
                {
                    return false;
                }
                BigInteger b = BigInteger.ModPow(a, t, R);
                if (b == 1 || b == R - 1)
                {
                    continue;
                }
                else
                    for (int j = 0; j < s - 1; j++)
                    {
                        b = (b * b) % R;
                        if (b == 1)
                        {
                            return false;
                        }
                        if (b == R - 1)
                        {
                            break;
                        }
                    }
                if (b == R - 1)
                {
                    continue;
                }
                return false;
            }
            return true;
        }
        public static BigInteger RandomBigInteger(BigInteger N)
        {
            byte[] bytes = N.ToByteArray();
            rand.NextBytes(bytes);
            bytes[bytes.Length - 1] &= 0x7F;
            BigInteger R = new BigInteger(bytes);
            while (R > N - 1)
            {
                R = R - N + 1;
            }
            if (R == 1)
            {
                R++;
            }
            if (R == 0)
            {
                R += 2;
            }
            return R;
        }
        public static BigInteger PrimeGenerator(int length)
        {
            BigInteger R = RandomBigInteger(BigInteger.Pow(2, length));
            if (R % 2 == 0)
                R++;
            int i = 2;
            while (true)
            {
                if (MillerRabinTest(R))
                    return R;
                else
                {
                    R = R + i;
                }
            }
        }

        public static BigInteger Invert(BigInteger a, BigInteger b)
        {
            BigInteger mod = b;
            if (a < 0)
                a += b;
            else if (a > b)
                a %= b;
            BigInteger q, r, x1, x2, y1, y2, x, y, inv, d;
            if (b == 0)
            {
                x = 1;
                y = 0;
                d = a;
                inv = x;
                return inv;
            }
            x1 = 0;
            x2 = 1;
            y1 = 1;
            y2 = 0;
            while (b > 0)
            {
                q = a / b;
                r = a - q * b;
                x = x2 - q * x1;
                y = y2 - q * y1;
                a = b;
                b = r;
                x2 = x1;
                x1 = x;
                y2 = y1;
                y1 = y;
            }
            x = x2;
            y = y2;
            d = a;
            inv = x;
            inv %= mod;
            if (inv < 0)
                inv += mod;
            return inv;
        }

        public static BigInteger[] Sieve(int n)
        {
            List<BigInteger> primes = new List<BigInteger>();
            bool[] isPrime = new bool[n];
            isPrime[0] = false;
            isPrime[1] = false;
            for (int i = 2; i < n; i++)
            {
                isPrime[i] = true;
            }
            for (int i = 0; i < n; i++)
            {
                if (isPrime[i])
                {
                    primes.Add(i);
                    for (int j = 2 * i; j < n; j += i)
                    {
                        isPrime[j] = false;
                    }
                }
            }
            return primes.ToArray();
        }
        public static BigInteger[] Sieve(int n, int m)
        {
            List<BigInteger> primes = new List<BigInteger>();
            bool[] isPrime = new bool[m];
            isPrime[0] = false;
            isPrime[1] = false;
            for (int i = 2; i < m; i++)
            {
                isPrime[i] = true;
            }
            for (int i = 2; i < n; i++)
            {
                if (isPrime[i])
                {
                    for (int j = 2 * i; j < m; j += i)
                    {
                        isPrime[j] = false;
                    }
                }
            }
            for (int i = n; i < m; i++)
            {
                if (isPrime[i])
                {
                    primes.Add(i);
                    for (int j = 2 * i; j < m; j += i)
                    {
                        isPrime[j] = false;
                    }
                }
            }
            return primes.ToArray();
        }

        public static BigInteger TonelliShanks(BigInteger n, BigInteger p)
        {
            //проверка числа
            if (LegendreSymbol(n, p) != 1)
            {
                return new BigInteger();
            }

            BigInteger q = p - 1;
            int s = 0;
            // Разложение числа 
            while ((q % 2) == 0)
            {
                q /= 2;
                s++;
            }
            if (s == 1)
                return BigInteger.ModPow(n, (p + 1) / 4, p);

            BigInteger z;
            Random random = new Random();
            // НАхождение квадратичного невычета z
            do
            {
                z = RandomBigInteger(p);
            }
            while (LegendreSymbol(z, p) != -1);
            // задание начальных значений 
            BigInteger c = BigInteger.ModPow(z, q, p);
            BigInteger result = BigInteger.ModPow(n, (q + 1) / 2, p);
            BigInteger t = BigInteger.ModPow(n, q, p);
            int m = s;
            while (t != 1)
            {
                BigInteger tt = t;
                int i = 1;
                for (; i < m; i++)
                {
                    tt = (tt * tt) % p;
                    if (tt == 1)
                        break;
                }
                BigInteger b = BigInteger.ModPow(c, BigInteger.Pow(2, m - i - 1), p);
                BigInteger b2 = (b * b) % p;
                result = (result * b) % p;
                t = (t * b2) % p;
                c = b2;
                m = i;
            }
            return result;
        }
        public static int LegendreSymbol(BigInteger a, BigInteger p)
        {
            if (a < 0)
                return -1;
            if (a % p == 0)
                return 0;
            if (a == 1)
                return 1;
            if (a % 2 == 0)
            {
                BigInteger temp = (p * p - 1) / 8;
                if (temp % 2 == 0)
                    return LegendreSymbol(a / 2, p);
                else
                    return LegendreSymbol(a / 2, p) * (-1);
            }
            else
            {
                BigInteger temp = (a - 1) * (p - 1) / 4;
                if (temp % 2 == 0)
                    return LegendreSymbol(p % a, a);
                else
                    return LegendreSymbol(p % a, a) * (-1);
            }
        }

        public static BigInteger GenerateRSA(int length1, int length2)
        {
            BigInteger p = PrimeGenerator(length1);
            BigInteger q = PrimeGenerator(length2);
            return p * q;
            //return BigInteger.Parse("41024879495073626381514254508702829516189620181801");
        }
    }
}
