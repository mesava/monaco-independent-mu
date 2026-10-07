# Независимая абсолютная нормировка Monte Carlo

Добавлен отдельный слой `AbsoluteMcCalibration`.

Он связывает:

1. **измеренную** reference dose/MU;
2. MC reference calculation в той же геометрии;
3. patient MC dose/history.

Если:

```text
D_meas_ref_per_MU = измеренная Gy/MU
D_MC_ref_per_history = MC Gy/history
```

то эквивалентное число source histories на MU:

```text
H_per_MU =
    D_meas_ref_per_MU / D_MC_ref_per_history
```

и для пациента:

```text
D_patient_per_MU =
    D_MC_patient_per_history * H_per_MU
```

Для конкретного beam:

```text
D_beam =
    D_patient_per_MU * BeamMeterset
```

## Критическое ограничение

Reference MC и patient MC должны иметь **одинаковую source-normalization
semantics**. Нельзя смешивать, например, разные определения incident history,
разные phase-space normalization или скрыто перенормированные outputs.

Поэтому в итоговом отчёте должны сохраняться:

- MC engine/version;
- source model hash;
- reference run hash;
- histories;
- scoring normalization;
- commissioning dataset hash.

## Почему Monaco reference сюда не входит

Значения 0.995/0.994/1.000/1.001 Gy за 100 MU из Monaco являются TPS reference
и используются только для последующего сравнения.

Если применить их вместо измеренной dose/MU для построения
`H_per_MU`, абсолютная проверка перестанет быть независимой.
