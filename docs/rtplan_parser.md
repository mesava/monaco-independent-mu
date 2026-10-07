# RTPLAN / VMAT parser

## Назначение

M2 преобразует DICOM RT Plan в явную модель delivery, не используя Monaco API
и не используя рассчитанную Monaco dose.

## BeamMeterset

MU не читается из `BeamSequence`. Он извлекается из:

```text
FractionGroupSequence
  -> ReferencedBeamSequence
      -> ReferencedBeamNumber
      -> BeamMeterset
```

и связывается с `BeamSequence` по `BeamNumber`.

Это принципиально важно для корректной DICOM RT semantics.

## Control points

Для каждого beam сохраняются:

- CumulativeMetersetWeight;
- gantry angle / rotation direction;
- collimator angle / rotation direction;
- patient support angle;
- tabletop positions;
- isocenter;
- NominalBeamEnergy;
- DoseRateSet;
- SSD;
- все BeamLimitingDevicePositionSequence.

Пропущенные в последующих CP атрибуты наследуют последнее явно заданное
состояние. Это относится и к jaw/MLC state.

## Beam limiting devices

Геометрия устройства берётся из `BeamLimitingDeviceSequence`:

- RTBeamLimitingDeviceType;
- NumberOfLeafJawPairs;
- LeafPositionBoundaries;
- SourceToBeamLimitingDeviceDistance.

Количество `LeafJawPositions` в каждом CP проверяется как
`2 * NumberOfLeafJawPairs`.

Для Agility ожидаем, что реальные данные определят тип MLC и число пар из
DICOM, поэтому код не hard-code'ит 80 пар / 160 координат.

## VMAT segment MU

Интервал delivery определяется соседними control points:

```text
delta_CMW = CMW[i+1] - CMW[i]

delta_MU =
    BeamMeterset * delta_CMW / FinalCumulativeMetersetWeight
```

Control points являются **границами динамического интервала**, а не статическими
полями, дозу которых следует просто суммировать. На этапе Monte Carlo геометрия
внутри сегмента должна интерполироваться/семплироваться по meterset weight.

## Safety checks

Parser останавливается при:

- нескольких Fraction Groups без явного выбора;
- отсутствующем BeamMeterset;
- неуникальном BeamNumber;
- undefined beam-limiting device;
- неправильном числе leaf/jaw positions;
- отсутствии начального состояния device;
- non-monotonic CMW;
- несовпадении последнего CMW и FinalCumulativeMetersetWeight.
