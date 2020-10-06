using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Numerics;

namespace ECM2018.core
{
    interface IEllipticCurve
    {
        BigInteger a { get;   }
        BigInteger b { get;   }
        BigInteger n { get;   }
    }
}
