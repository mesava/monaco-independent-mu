# Независимый расчёт MU для Monaco

Исследовательский проект независимой проверки мониторных единиц (MU) и трёхмерного распределения дозы для планов лучевой терапии, экспортированных из **Elekta Monaco 6.1.4**, с использованием **Monte Carlo**, DICOM RT и пациентской КТ.

> ⚠️ **Статус:** исследовательская разработка. Проект не предназначен для клинического применения без полного commissioning, валидации, контроля версий и локального утверждения методики.

## Цель

Создать независимый расчётный контур для Elekta Versa HD / Agility, который:

- читает RT Plan, RT Dose, RT Structure Set и CT;
- восстанавливает геометрию каждого beam/control point;
- учитывает динамику MLC Agility и VMAT;
- строит пациентскую модель по CT с учётом гетерогенностей;
- выполняет независимый Monte Carlo расчёт абсолютной дозы;
- выполняет **per-beam equivalent-MU consistency check** по независимой dose-per-MU;
- использует несколько малых ROI для устойчивого per-beam MU-check;
- выполняет независимое 3D сравнение MC vs Monaco;
- формирует автоматические QA-флаги и отчёт.

## Зафиксированные критерии согласия

Для MU-check:

- простые планы (open fields, tangents): **±3%**;
- IMRT/VMAT: **±5%**;
- превышение установленного порога автоматически отмечается флагом;
- итоговая градация: **GREEN / YELLOW / RED**.

Точные правила GREEN/YELLOW/RED будут окончательно утверждены после валидационной серии и не должны рассматриваться как универсальные клинические допуски до её завершения.

## Текущая архитектура

```text
Monaco / DICOM
      │
      ▼
DICOM parser
  ├── RT Plan
  ├── RT Dose
  ├── RT Struct
  └── CT
      │
      ▼
Patient model
  HU → RED → mass density → material
      │
      ├──────────────┐
      ▼              ▼
Beam model       MLC Agility
Versa HD         control points
      │              │
      └──────┬───────┘
             ▼
       Monte Carlo
      EGSnrc / GPU
             │
      ┌──────┴───────┐
      ▼              ▼
per-beam MU      3D dose
multi-ROI        gamma/DVH
      │              │
      └──────┬───────┘
             ▼
      QA / отчёт
```

## Текущее состояние

Проект уже вышел за рамки первоначального patient-model prototype.

| Milestone | Состояние |
|---|---|
| **M1 — Patient model** | HU→RED→ρ→MaterialMix, PEGSless, egsphant, sensitivity framework реализованы |
| **M2 — RTPLAN/VMAT** | BeamMeterset, CP inheritance, jaws/MLC, CMW и delivery segments реализованы |
| **M3 — Versa HD/Agility** | commissioning data model, TPS reference, machine/energy resolution, IEC bank identity и research rounded-tip transform реализованы; physical head model ещё валидируется |
| **M4 — Monte Carlo** | transport abstraction, SYNCMLCE/SYNCJAWS/source21 layers, 3ddose и absolute-normalization infrastructure реализованы; полноценный commissioned transport model ещё не закрыт |
| **M5 — Equivalent MU** | multi-ROI equivalent-MU check, patient-coordinate spherical ROI и ±3%/±5% infrastructure реализованы |
| **M6 — 3D comparison** | RTDOSE/MC alignment, ΔD, DVH и explicit gamma wrapper реализованы |

### Основные открытые физические задачи

До первого полноценного независимого расчёта клинического плана должны быть
закрыты четыре принципиальных узла:

1. **Agility rounded-tip geometry.** Реализованы analytic tangent transform,
   IEC bank identity и research-only SYNCMLCE mapper; остаётся end-to-end
   validation на реальной SYNCMLCE geometry с eccentric tip/leaf-bank tilt.
2. **Source 21 orientation.** Реализован HFS IEC/DICOM → DOSXYZnrc transform
   Zhan/Jiang/Osei с non-coplanar regression tests; остаётся end-to-end
   transport validation в DOSXYZnrc/BEAMnrc.
3. **Независимая absolute calibration.** Monaco reference doses за 100 MU
   хранятся только как TPS benchmark; MC нормируется по измеренной dose/MU.
4. **Transport commissioning.** Open fields → MLC stress tests → IMRT → VMAT
   должны пройти сравнение с измерениями до клинического использования.

### Реальный validation case 001

Получен полный Monaco DICOM export и выполнен технический preflight без
коммита исходных пациентских данных в Git.

Подтверждено:

- 187 CT slices, HFS, 120 kV;
- 8-field 6-MV dynamic-MLC IMRT;
- Agility MLCX: 80 pairs × 5 mm, span 400 mm;
- plan RTDOSE в Gy;
- RTSTRUCT без external/BODY/patient ROI.

Последний пункт выявил реальный workflow gap. Для него добавлен только
явно выбираемый research CT-derived external-mask path, seeded by treatment
isocenter. Он не является клиническим default и должен быть независимо
проверен до использования.

Подробности: docs/validation_case_001_dicom.md.

### Принцип независимости

```text
Измеренная commissioning dose/MU
             │
             ▼
Independent MC ───────────────► Gy/MU
             │
             ├──► equivalent MU per beam
             │
             └──► 3D dose / gamma / DVH
                         ▲
                         │
                  Monaco RTDOSE
                  только reference
```

Monaco dose **никогда не используется для подгонки или абсолютной нормировки
independent MC**.

## План разработки

1. **M1 — Patient model**
   - DICOM CT geometry;
   - HU → RED;
   - mass density;
   - patient material assignment;
   - экспорт voxel phantom.

2. **M2 — RTPLAN / VMAT parser**
   - beam geometry;
   - MU;
   - control points;
   - gantry/collimator/couch;
   - jaws;
   - MLC leaf positions;
   - CumulativeMetersetWeight.

3. **M3 — Versa HD / Agility beam model**
   - commissioning PDD/profiles/output factors;
   - абсолютная калибровка;
   - Agility geometry;
   - transmission/leakage;
   - leaf offset/gap;
   - FFF/flattened beams.

4. **M4 — Monte Carlo engine**
   - EGSnrc-compatible backend;
   - VMAT integration;
   - GPU acceleration;
   - статистическая неопределённость.

5. **M5 — Independent MU**
   - per-beam MU;
   - multi-ROI;
   - ROI radius: **0.25 cm**;
   - mean dose in ROI;
   - GREEN/YELLOW/RED.

6. **M6 — 3D comparison**
   - абсолютная MC dose grid;
   - приведение к сетке RT Dose;
   - 3D gamma;
   - DVH comparison;
   - voxel-wise dose difference.

## Язык

Основное вычислительное ядро и orchestration: **Python**.

Причины:

- pydicom + NumPy/SciPy;
- удобная работа с 3D массивами;
- исследовательская обработка commissioning data;
- запуск внешнего MC backend;
- GPU tooling;
- тестируемость и быстрые итерации.

C# может быть добавлен позднее как Windows/clinical frontend или слой интеграции с Monaco API.

## Структура репозитория

```text
src/indep_mu/
├── dicom/
├── patient_model/
├── beam_model/
├── montecarlo/
├── mu_check/
└── comparison3d/

config/
├── ct/
├── materials/
└── machine/

docs/
tests/
data/
```

## Данные пациентов

**Пациентские DICOM не должны коммититься в репозиторий.**

В Git хранятся:

- исходный код;
- обезличенные конфигурации;
- machine/commissioning data без персональных данных;
- тестовые синтетические данные;
- агрегированные результаты валидации.

В `.gitignore` блокируются DICOM, архивы, MC dose files и каталоги пациентских данных.

## Требование независимости

Monaco используется как источник:

- DICOM;
- фактических MU;
- исходного RT Dose для сравнения;
- scanner-specific CT calibration;
- machine commissioning/configuration data.

Независимый движок **не должен использовать Monaco dose engine для вычисления независимой дозы**.

## Валидация

До любого клинического применения необходимы как минимум:

- unit tests геометрии и DICOM;
- тесты HU/RED/density conversion;
- водный/гомогенный фантом;
- heterogeneous phantom;
- static fields;
- MLC stress tests;
- IMRT;
- VMAT;
- сравнение с измерениями;
- серия клинических планов;
- анализ статистической и систематической неопределённости.

---

Проект ведётся последовательно: каждый физический слой должен быть валидирован до перехода к следующему.
