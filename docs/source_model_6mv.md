# 6 MV source-model commissioning strategy

## Почему сейчас фокус только на 6 MV

Первый реальный validation case использует только **6 MV STANDARD**. Поэтому
для первого end-to-end независимого расчёта нет необходимости одновременно
закрывать 6 FFF, 10 MV и 10 FFF.

Это позволяет выстроить commissioning последовательно:

1. 6 MV static open fields;
2. 6 MV Agility MLC tests;
3. 6 MV dynamic-MLC IMRT validation case 001;
4. только после этого переносить framework на остальные energy models.

## Literature priors — не финальная модель

Paschal et al. моделировали 6-MV Versa HD в BEAMnrc с Source 19 и использовали
incident electron energy около **6.4 MeV**.

Gholampourkashi et al. для 6-MV Elekta Infinity / Agility получили tuned
source parameters:

- electron energy: **6.6 ± 0.1 MeV**;
- cross-line beam width: **2.1 ± 0.1 mm**;
- in-line beam width: **1.0 ± 0.1 mm**;
- angular divergence: **1.35 ± 0.20°**.

В той же работе source energy подбиралась по PDD, а source widths — по
измеренным penumbrae/profiles.

Эти числа являются только стартовой областью поиска. Они не являются
commissioning data нашей машины.

## Local commissioning rule

Параметры должны определяться по **локальным измерениям**, а не по Monaco dose:

- electron energy → PDD;
- cross-line source width → cross-line penumbra/profile;
- in-line source width → in-line penumbra/profile;
- angular divergence → profile/field-size dependence;
- absolute scale → measured reference dose/MU.

Monaco допускается только как independent comparison target после того, как
source/head model уже определён по измерениям.

## Agility parameters

Gholampourkashi также получил для своей модели:

- LBROT ≈ 9 mrad;
- leaf density ≈ 18.5 g/cm³ after tuning;
- average measured leaf transmission ≈ 4.3%;
- MC leaf transmission ≈ 4.1%;
- nominal interleaf air gap ≈ 0.089 mm.

Эти параметры полезны как sensitivity priors, но density/composition/gap не
должны переноситься в нашу модель автоматически.

## First 6-MV validation ladder

Перед patient case 001:

1. 10×10 reference field;
2. PDD 5×5 / 10×10 / larger fields;
3. cross-line and in-line profiles;
4. output factors;
5. MLC transmission;
6. alternating open/closed leaves / FOURL-like pattern;
7. dynamic/sweeping-gap delivery;
8. только затем 8-field DMLC IMRT patient plan.

## Sources

- Paschal HMP et al. J Appl Clin Med Phys. 2022;23:e13715.
  DOI: 10.1002/acm2.13715.
- Gholampourkashi S et al. J Appl Clin Med Phys. 2019;20:55–67.
  DOI: 10.1002/acm2.12485.
