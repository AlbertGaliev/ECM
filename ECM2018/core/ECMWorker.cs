using ECM2018.core.weierstrass;
using ECM2018.core.montgomery;
using ECM2018.core.edwards;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;
using ECM2018.core.algorithms;
using System.Threading;

namespace ECM2018.core
{
    class ECMWorker
    {
        object locker;
        bool stopLoop;
        bool stopFactoring;
        public List<BigInteger> primeFactors = new List<BigInteger>();
        public List<BigInteger> compositeFactors = new List<BigInteger>();
        static Random rand = new Random();
        ECMParams param;
        public ECMWorker(ECMParams param)
        {
            this.param = param;
        }

        public BigInteger UseMethod(BigInteger n)
        {
            //if (n< 1000000)
            //{
            //    return ECMWeierstrass.GetFactorBothStages(n, param.B1, param.B2);
            //}
            switch (param.curveType)
            {
                case CurveType.WeierstrassCurve:
                    switch (param.coordinatesType)
                    {
                        case CoordinatesType.Affine:
                            return ECMWeierstrass.GetFactorAffine(n, param.B1, param.B2);
                        case CoordinatesType.Projective:
                            return ECMWeierstrass.GetFactorProjective(n, param.B1, param.B2);
                        case CoordinatesType.Jacobian:
                            return ECMWeierstrass.GetFactorJacobian(n, param.B1, param.B2);
                    }
                    break;
                case CurveType.MontgomeryCurve:
                    switch (param.coordinatesType)
                    {
                        case CoordinatesType.Projective:
                            return ECMMontgomery.GetFactorProjective(n, param.B1, param.B2, param.torsionType);
                    }
                    break;
                case CurveType.TwistedEdwardsCurve:
                    switch (param.coordinatesType)
                    {
                        case CoordinatesType.Projective:
                            return ECMEdwards.GetFactorProjective(n, param.B1, param.B2, param.torsionType);
                        case CoordinatesType.Inverted:
                            return ECMEdwards.GetFactorInverted(n, param.B1, param.B2, param.torsionType);
                        case CoordinatesType.Extended:
                            return ECMEdwards.GetFactorExtended(n, param.B1, param.B2, param.torsionType);
                            //case CoordinatesType.Completed:
                            //    return ECMEdwards.GetFactorProjective(n, param.B1, param.B2, param.torsionType);
                    }
                    break;
            }
            return 0;
        }
        public BigInteger UseMethod(BigInteger n, int b1, int b2)
        {
            if (n < 1000000)
            {
                return ECMWeierstrass.GetFactorAffine(n, b1, b2);
            }
            switch (param.curveType)
            {
                case CurveType.WeierstrassCurve:
                    switch (param.coordinatesType)
                    {
                        case CoordinatesType.Affine:
                            return ECMWeierstrass.GetFactorAffine(n, b1, b2);
                        case CoordinatesType.Projective:
                            return ECMWeierstrass.GetFactorProjective(n, b1, b2);
                        case CoordinatesType.Jacobian:
                            return ECMWeierstrass.GetFactorJacobian(n, b1, b2);
                    }
                    break;
                case CurveType.MontgomeryCurve:
                    switch (param.coordinatesType)
                    {
                        case CoordinatesType.Projective:
                            return ECMMontgomery.GetFactorProjective(n, b1, b2, param.torsionType);
                    }
                    break;
                case CurveType.TwistedEdwardsCurve:
                    switch (param.coordinatesType)
                    {
                        case CoordinatesType.Projective:
                            return ECMEdwards.GetFactorProjective(n, b1, b2, param.torsionType);
                        case CoordinatesType.Inverted:
                            return ECMEdwards.GetFactorInverted(n, b1, b2, param.torsionType);
                        case CoordinatesType.Extended:
                            return ECMEdwards.GetFactorExtended(n, b1, b2, param.torsionType);
                            //case CoordinatesType.Completed:
                            //    return ECMEdwards.GetFactorProjectiveLoop(n, param.torsionType);
                    }
                    break;
            }
            return 0;
        }

        public BigInteger GetFactor(BigInteger n)
        {
            ECM.stop = false;
            locker = new object();
            BigInteger currentFactor = 0;
            if (!param.autoApproximation)
            {
                Parallel.For(0, param.curves, (number, state) =>
                {
                    var result = UseMethod(n);
                    if (result != 0)
                        lock (locker)
                        {
                            currentFactor = result;
                            state.Stop();
                        }
                });
            }
            else
            {
                stopLoop = false;
                Parallel.For(0, param.curves, (number, state) =>
                {
                    int B1 = 10;
                    int B2 = 20;
                    int iterations = 1;
                    while (!stopLoop)
                    {
                        for (int k = 0; k < iterations; k++)
                        {
                            var result = UseMethod(n, B1, B2);
                            if (stopLoop)
                                break;
                            if (result > 1)
                                lock (locker)
                                {
                                    currentFactor = result;
                                    stopLoop = true;
                                    state.Stop();
                                    break;
                                }
                        }
                        B1 *= 3;
                        B2 *= 3;
                        iterations *= 2;
                    }
                });
            }
            return currentFactor;
        }
        public BigInteger[] GetAllFactors(BigInteger n)
        {
            stopFactoring = false;
            compositeFactors.Clear();
            primeFactors.Clear();
            stopLoop = false;
            if (Algorithms.MillerRabinTest(n))
            {
                primeFactors.Add(n);
                return primeFactors.ToArray();
            }
            else
                compositeFactors.Add(n);
            while (compositeFactors.Count > 0)
            {
                if (stopFactoring)
                    break;
                var number = compositeFactors[0];
                compositeFactors.Remove(number);
                var factor1 = GetFactor(number);
                if (factor1 == 0)
                {
                    compositeFactors.Add(number);
                    continue;
                }
                var factor2 = number / factor1;
                if (Algorithms.MillerRabinTest(factor1))
                    primeFactors.Add(factor1);
                else
                    compositeFactors.Add(factor1);
                if (Algorithms.MillerRabinTest(factor2))
                    primeFactors.Add(factor2);
                else
                    compositeFactors.Add(factor2);
            }
            primeFactors.Sort();
            return primeFactors.ToArray();
        }
        public void StopFactoring()
        {
            ECM.stop = true;
            stopLoop = true;
            stopFactoring = true;
        }


        public class ECMParams
        {
            public CurveType curveType { get; set; }
            public CoordinatesType coordinatesType { get; set; }
            public TorsionType torsionType { get; set; }
            public bool autoApproximation { get; set; }
            public int B1 { get; set; }
            public int B2 { get; set; }
            public int curves { get; set; }
            public ECMParams()
            {

            }
            public ECMParams(CurveType curve, CoordinatesType coordinates, TorsionType torsion,
                bool autoApprox, int curves, int B1, int B2)
            {
                this.curveType = curve;
                this.coordinatesType = coordinates;
                this.torsionType = torsion;
                this.autoApproximation = autoApprox;
                this.curves = curves;
                this.B1 = B1;
                this.B2 = B2;
            }
        }

    }
}
