#pragma once

// SPICE-compatible setup
class SPICECompatible{
public:
    // Default constructor that zero-initializes or sets desired defaults
    SPICECompatible() = default;
    
    /* defines for CKTmode */
    enum SPICEmode : unsigned long {
        /* old 'mode' parameters */
        MODE               = 0x3,
        MODETRAN           = 0x1,
        MODEAC             = 0x2,

        /* for noise analysis */
        MODEACNOISE        = 0x8,

        /* old 'modedc' parameters */
        MODEDC            = 0x70,
        MODEDCOP          = 0x10,
        MODETRANOP        = 0x20,
        MODEDCTRANCURVE   = 0x40,

        /* old 'initf' parameters */
        INITF           = 0x3f00,
        MODEINITFLOAT    = 0x100,
        MODEINITJCT      = 0x200,
        MODEINITFIX      = 0x400,
        MODEINITSMSIG    = 0x800,
        MODEINITTRAN    = 0x1000,
        MODEINITPRED    = 0x2000,  

        /* old 'nosolv' paramater */
        MODEUIC          = 0x10000l
    };
    
    inline void setMode(unsigned long spiceCompatibleMode)
    {
        CKTmode = spiceCompatibleMode;
    }
    inline unsigned long getMode() const
    {
        return CKTmode;
    }

    void setFlagsOP();
    void setFlagsTranOP();
    void setFlagsTR();
    void setFlagsDC();
    void setFlagsSmallSig();
    void setFlagsAC();
    void updateStateMachine(bool converged, bool &NISHOULDREORDER, int iterno);

private:
    unsigned long CKTmode;
};

// .op
void SPICECompatible::setFlagsOP(){
    setMode((CKTmode & MODEUIC) | MODEDCOP | MODEINITJCT);
}
// .tran
void SPICECompatible::setFlagsTranOP(){
    setMode((CKTmode & MODEUIC) | MODETRANOP | MODEINITJCT);
}
void SPICECompatible::setFlagsTR(){
    setMode((CKTmode & MODEUIC) | MODETRAN | MODEINITTRAN);
}
// .dc
void SPICECompatible::setFlagsDC(){
    setMode((CKTmode & MODEUIC) | MODEDCTRANCURVE | MODEINITJCT);
}
// .ac
void SPICECompatible::setFlagsSmallSig(){
    setMode((CKTmode & MODEUIC) | MODEDCOP | MODEINITSMSIG);
}
void SPICECompatible::setFlagsAC(){
    setMode((CKTmode & MODEUIC) | MODEAC);
}

void SPICECompatible::updateStateMachine(bool converged, bool &NISHOULDREORDER, int iterno)
{
    if (CKTmode & MODEAC)
    {
    }
    else if (CKTmode & MODEINITFLOAT)
    {
        /*if ((CKTmode & MODEDC) && ckt->CKThadNodeset)
        {
            if (ipass)
                ckt->CKTnoncon = ipass;
            ipass = 0;
        }*/
    }
    else if (CKTmode & MODEINITJCT)
    {
        CKTmode = (CKTmode & ~INITF) | MODEINITFIX;
        /*ckt->CKTniState |= NISHOULDREORDER;*/
        NISHOULDREORDER = true;
    }
    else if (CKTmode & MODEINITFIX)
    {
        if (converged)
            CKTmode = (CKTmode & ~INITF) | MODEINITFLOAT;
    }
    else if (CKTmode & MODEINITSMSIG)
    {
        CKTmode = (CKTmode & ~INITF) | MODEINITFLOAT;
    }
    else if (CKTmode & MODEINITTRAN)
    {
        if (iterno <= 1)
            NISHOULDREORDER = true;
            /*ckt->CKTniState |= NISHOULDREORDER;*/
        CKTmode = (CKTmode & ~INITF) | MODEINITFLOAT;
    }
    else if (CKTmode & MODEINITPRED)
    {
        CKTmode = (CKTmode & ~INITF) | MODEINITFLOAT;
    }

}
