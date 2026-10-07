# Research geometry: projected Agility leaf edge to SYNCMLCE cylinder origin

## Why simple ZMIN/SAD scaling is wrong for a rounded tip

BEAMnrc SYNCMLCE uses different meanings for NEG/POS.

For ENDTYPE=1 they are coordinates of the opening at ZMIN. For ENDTYPE=0
they are transverse coordinates of the origin of the cylinder that defines the
rounded leaf end.

The SYNCMLCE source code evaluates a rounded tip through the distance from the
particle point to a circle with centre (leaf_origin, CIL) and radius R. Thus
LEAFB/LEAFA are cylinder-centre coordinates and CIL is the z coordinate of the
centre of curvature.

## Idealized tangent mapping

Place the source at (0,0), the isocenter plane at z=SAD and let x_iso be the
projected leaf edge at isocenter.

The source ray is

    x = m z
    m = x_iso / SAD

For a cylinder centred at (c, CIL), tangency requires

    abs(c - m*CIL) / sqrt(1 + m^2) = R

Therefore

    negative opening side:
    c = m*CIL - R*sqrt(1+m^2)

    positive opening side:
    c = m*CIL + R*sqrt(1+m^2)

The project now contains RoundedLeafTipTangentGeometry, which implements both
directions of this transform and tests exact tangency and round-trip behaviour.

## Public Agility starting values

The research helper uses only geometric values reported in public literature:

    R   = 17.0 cm
    CIL = 34.93 cm
    SAD = 100 cm by default

Hernandez et al. reports 34.93 cm as source-to-collimator distance measured at
the centre of curvature and a rounded-tip radius of 17 cm. These values are
consistent with other published Agility Monte Carlo models.

## Why this is not yet a production DICOM-to-SYNCMLCE mapper

The equation above describes an ideal central cylindrical tip before the rest
of SYNCMLCE geometry is applied.

SYNCMLCE additionally:

1. creates leaves with divergence/focusing according to leaf-centre spacing;
2. translates and rotates each leaf relative to the central leaf;
3. applies LBROT to the complete leaf bank;
4. applies step/tongue-and-groove geometry.

The project therefore does not yet feed DICOM Leaf/Jaw Positions through this
formula directly into a clinical Agility sequence.

The DICOM bank order is no longer left ambiguous: DICOM uses IEC element
order 101..1N followed by 201..2N, and IEC 61217 defines the 1-side
(X1/Y1) as the negative-axis edge and the 2-side (X2/Y2) as the
positive-axis edge. The code therefore has a standards-based bank identity
helper; this still needs vendor-export confirmation on the first real Monaco
RTPLAN.

The research mapper also accepts an explicit projected-edge shift. This is
intended for the field-centre translation caused by leaf-bank tilt. The code
does not infer the sign: the sign must be fixed by the validated machine-axis
convention. For a 90 mm leaf thickness and LBROT=9 mrad the published
magnitude is about 0.405 mm.

The remaining items to verify are:

- exact Monaco/DICOM meaning of the nominal edge for eccentric Agility tips;
- sign and application point of the LBROT field-centre shift;
- interaction of the tangent transform with per-leaf focusing;
- tilted leaf-side/interleaf geometry;
- agreement of the actual SYNCMLCE surface with projected field edges.

## Next validation step

For a static central leaf pair:

1. compare analytic source-ray tangency with the actual SYNCMLCE surface;
2. recover the projected field edge at isocenter;
3. introduce LBROT=9 mrad and quantify the induced shift;
4. compare with published/measurement-derived translation;
5. only then extend the mapper to all 80 leaf pairs and dynamic sequences.

## Sources

- NRC EGSnrc/BEAMnrc SYNCMLCE input documentation and source code.
- Hernandez V et al. Medical Physics 2022. DOI 10.1002/mp.16016.
- Gholampourkashi S et al. JACMP 2019. DOI 10.1002/acm2.12485.
- Ohira S et al. In Vivo 2020. DOI 10.21873/invivo.12050.

The research mapper is accepted by build_syncmlce_sequence only when selected
explicitly by the caller. No clinical/default configuration selects it, and
head-geometry readiness still blocks production use until transport validation
is complete.
