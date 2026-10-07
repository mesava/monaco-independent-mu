# Публичная геометрия Agility: исследовательский baseline

Для продвижения независимой модели создан отдельный
`head_geometry.literature_candidate.yaml`.

Он **не является клинически утверждённой геометрией Versa HD**. Цель файла —
отделить публично подтверждаемые размеры от неизвестных параметров и сделать
все допущения видимыми.

## Что подтверждается открытой литературой

Hernandez et al. (Medical Physics, 2022) приводит для Agility:

- 160 leaves / 80 pairs;
- projected leaf width 5 mm at isocenter;
- leaf thickness 9.0 cm;
- source-to-collimator distance 34.93 cm, определённую в центре кривизны;
- rounded-tip radius 17 cm;
- eccentric rounded tip.

В той же работе центр кривизны описан как расположенный:

- 37.5 mm от верхней поверхности leaf;
- 52.5 mm от нижней поверхности leaf.

Отсюда для исследовательской геометрии:

```text
z_center = 34.93 cm
zmin = 34.93 - 3.75 = 31.18 cm
zmax = 34.93 + 5.25 = 40.18 cm
```

Gholampourkashi et al. получил наилучшее согласие с измерениями для
`LBROT ≈ 9 mrad`; Hernandez et al. цитирует эквивалентное значение около
0.515°.

Ohira et al. независимо использовал:

- rounded-tip radius 17 cm;
- leaf thickness 9 cm;
- tungsten alloy density 18.0 g/cm³;
- composition 95% W, 3.75% Ni, 1.25% Fe;
- interleaf air gap 0.009 cm.

Эти material values пока не объявляются единственными правильными: другая
публикация показала чувствительность transmission к эффективной плотности и
использовала tuning против измерений.

## Что НЕ считается установленным

Пока остаются заблокированными:

- end-to-end validation nominal DICOM leaf position → фактическая
  SYNCMLCE rounded-tip surface; analytic tangent transform уже реализован;
- знак/применение компенсации field-centre shift от LBROT;
- полная геометрия Y-jaws, включая нижнюю поверхность;
- material density, окончательно выбранная по transmission measurements;
- остальные proprietary head components для полного Versa HD.

DICOM/IEC bank identity теперь не считается неизвестной: DICOM хранит позиции
в IEC element order 101..1N, 201..2N, а IEC 61217 определяет side 1
(X1/Y1) как negative-axis side и side 2 (X2/Y2) как positive-axis side.
На первом реальном Monaco RTPLAN это всё равно будет подтверждено
sanity-check'ом экспортированных позиций.

## Почему public model всё равно полезен

Он задаёт независимый starting point для последующей commissioning-driven
оптимизации. Параметры не должны подгоняться к Monaco dose. Они будут
проверяться против измеренных:

- PDD;
- profiles;
- output factors;
- MLC transmission;
- FOURL / alternating leaf patterns;
- sweeping-gap / dynamic tests.

## Источники

- Hernandez V et al. *Challenges in modeling the Agility multileaf collimator
  in treatment planning systems and current needs for improvement.* Med Phys.
  2022. DOI: 10.1002/mp.16016.
- Gholampourkashi S et al. *Monte Carlo and analytic modeling of an Elekta
  Infinity linac with Agility MLC.* J Appl Clin Med Phys. 2019.
  DOI: 10.1002/acm2.12485.
- Ohira S et al. *Monte Carlo Modeling of the Agility MLC for IMRT and VMAT
  Calculations.* In Vivo. 2020. DOI: 10.21873/invivo.12050.
- Paschal HMP et al. *Monte Carlo modeling of the Elekta Versa HD and patient
  dose calculation with EGSnrc/BEAMnrc.* J Appl Clin Med Phys. 2022.
  DOI: 10.1002/acm2.13715.


## Проверка против реального Monaco RTPLAN

Validation case 001 экспортирует для MLCX:

    SourceToBeamLimitingDeviceDistance = 349 mm

Это находится всего на 0.3 mm от опубликованного расстояния до центра
кривизны Agility rounded tip:

    CIL = 349.3 mm

Такое согласие является полезной независимой cross-check: DICOM reference
plane MLC практически совпадает с опубликованным centre-of-curvature plane.

Однако DICOM определяет (300A,00BA) как расстояние от radiation source до
beam-limiting device в целом. Этот tag сам по себе не определяет upstream/downstream
physical surface и поэтому **не заменяет** zmin/zmax физической leaf geometry.

В validator добавлена проверка MLC source distance около 349.3 mm с
1-mm tolerance. Она проверяет delivery geometry, но не объявляет DICOM
distance физической толщиной leaf.
