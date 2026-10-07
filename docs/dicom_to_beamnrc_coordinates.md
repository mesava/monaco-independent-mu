# DICOM IEC coordinates → BEAMnrc synchronized components

## Почему этот слой отделён

DICOM и BEAMnrc описывают положение ограничителей **не одной и той же
величиной**.

DICOM `Leaf/Jaw Positions (300A,011C)` содержит 2N координат в IEC BEAM
LIMITING DEVICE axis в порядке:

```text
101, 102, ..., 1N, 201, 202, ..., 2N
```

Поэтому parser сохраняет bank 1 и bank 2 отдельно. Код намеренно не угадывает
bank identity по знаку координаты.

`Leaf Position Boundaries (300A,00BE)` являются механическими границами
листьев, спроецированными в изоцентрическую плоскость. DICOM также использует
проекцию leaf/jaw geometry к isocenter plane в соответствующей модели beam
limiting device.

## Простая source projection

Если coordinate подтверждён как isocenter-projected и физическая поверхность
сфокусирована в источник, то:

```text
x(z) = x_iso * z / SAD
```

В проекте это реализовано только через явный объект `IsocenterProjection`.
Никакого скрытого масштабирования в DICOM parser нет.

## SYNCMLCE

BEAMnrc различает два принципиально разных случая:

### ENDTYPE = 1 — focused divergent leaf end

`NEG/POS` — координаты **front opening at ZMIN**.

Для leaf end, сфокусированного в источник, допустимо прямое преобразование:

```text
opening_ZMIN = opening_iso * ZMIN / SAD
```

Для этого реализован `SourceFocusedSyncMlceMapper`.

### ENDTYPE = 0 — rounded/cylindrical leaf end

`NEG/POS` означают уже **origin of the cylinder defining the rounded leaf
end**, а не projected field edge.

Поэтому DICOM leaf position нельзя просто умножить на `ZMIN/SAD` и записать в
SYNCMLCE sequence.

Production Agility mapping для rounded tips будет включён только после
фиксации и валидации полной leaf-tip geometry:

- radius;
- CIL / cylinder-axis z;
- bank tilt / LBROT;
- relation between DICOM projected leaf edge and cylinder origin.

До этого проект технически не позволяет замаскировать простую проекцию под
Agility rounded-tip model.

## Bank mapping

BEAMnrc требует `NEG` и `POS`. DICOM хранит bank 1 и bank 2.

Какой DICOM bank соответствует отрицательной стороне, задаётся **явно**
(`negative_bank=1|2`) в head-model configuration. Автоматическое определение
по знаку запрещено, потому что динамическая MLC может пересекать центральную
ось.

## MUINDEX

Для каждого beam проект запускает независимый MC отдельно. Поэтому:

```text
MUINDEX = CumulativeMetersetWeight / FinalCumulativeMetersetWeight
```

и последовательность идёт 0 → 1 внутри beam.

Абсолютные `BeamMeterset` используются позже при Gy/MU normalization и
суммировании beam doses.

Если несколько соседних DICOM CP имеют одинаковый CMW, между ними доставлено
0 MU. Для synchronized sequence сохраняется последнее состояние этого
zero-weight transition, чтобы исключить деление на ноль при динамической
интерполяции.

## SYNCJAWS

Для source-focused jaws физические координаты front/back поверхностей
вычисляются независимо:

```text
X_front = X_iso * ZMIN / SAD
X_back  = X_iso * ZMAX / SAD
```

Sequence writer сохраняет формат BEAMnrc:

```text
NFIELDS
INDEX
ZMIN ZMAX XFP XBP XFN XBN
...
```

Геометрические `ZMIN/ZMAX` не берутся автоматически из DICOM
`SourceToBeamLimitingDeviceDistance`: они являются частью валидированной
physical head model.

## Источники

- DICOM PS3.3, RT Beams / Beam Limiting Device Position Macro.
- DICOM PS3.3 C.8.8.25.3, Leaf Position Boundaries.
- EGSnrc PIRS-509A, SYNCMLCE input format.
- EGSnrc `SYNCJAWS_cm.mortran`.
