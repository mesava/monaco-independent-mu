# Физические источники и зафиксированные правила

## Monaco CT → RED

Используется scanner-specific таблица, предоставленная для текущей модели Monaco:
`DICOM3.DRT120kV`.

Elekta указывает, что значения CT ниже/выше диапазона CT-to-ED таблицы получают
соответственно минимальный/максимальный RED таблицы. В независимом движке это
реализуется как явный режим endpoint clipping, а safety layer всё равно обязан
сообщать количество out-of-range voxels.

Официальный источник:

- Elekta, *Monaco Physics Training — CT-ED*:
  https://www.elekta.co.jp/services/software-download/monaco/166/Monaco_Physics_Training_web_CT-ED_20201202.pdf

## RED → mass density

Используется piecewise relation Monaco:

- RED <= 0: rho = 0;
- 0 < RED < 1:
  rho = (sqrt(0.99^2 + 4*0.01*RED) - 0.99) / (2*0.01);
- RED >= 1:
  rho = (RED - 0.15) / 0.85.

Официальный источник:

- Elekta, *Patient Model Review*:
  https://www.elekta.co.jp/services/software-download/unity/217/6-1_Patient_Model_Review5.51_01_Japan_Web.pdf

## Mass density → material

Для Patient и Phantom LUT Monaco интерполирует материал между соседними
плотностями LUT. Например, диапазон 1.08–1.85 g/cm3 представляет смесь
MuscleSkeletalIcrp ↔ BoneCorticalIcrp.

Для PPS LUT Monaco не смешивает материалы между соседними точками.

Источник: тот же *Patient Model Review*, раздел CT to Mass Density.

## Независимая QA-политика

Даже если Monaco-matched calculation использует endpoint clipping, независимый
preflight сохраняет out-of-range voxels как отдельный QA-сигнал. Никакое
clipping не должно происходить молча.


## Independent elemental compositions

Для EGSnrc reference backend базовые тканевые composition берутся независимо
от Monaco из NIST STAR:

- Dry Air:
  https://physics.nist.gov/PhysRefData/Star/Text/table2.html
- Adipose Tissue (ICRP):
  https://physics.nist.gov/cgi-bin/Star/compos.pl?matno=103
- Muscle, Skeletal (ICRP):
  https://physics.nist.gov/cgi-bin/Star/compos.pl?matno=201
- Bone, Cortical (ICRP):
  https://physics.nist.gov/cgi-bin/Star/compos.pl?matno=120

EGSnrc `material.dat` независимо подтверждает standard dry-air and ICRP
cortical-bone definitions:
https://github.com/nrc-cnrc/EGSnrc/blob/master/HEN_HOUSE/pegs4/data/material.dat
