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

private:
    unsigned long CKTmode;
};

// .op
void SPICECompatible::setFlagsOP(){


}
// .tran
void SPICECompatible::setFlagsTranOP(){
    setMode((CKTmode & MODEUIC) | MODETRANOP | MODEINITJCT);
}

