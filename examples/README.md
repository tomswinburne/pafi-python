# PAFI Examples

## Systems
We have two point defects relaxed with an EAM potential, 
with another potential for testing thermodynamic integration: 
- `systems/EAM-EAM-SIA-Fe` : <br>
    dumbell SIA with Marinica07 EAM Fe (base) and Mendelev EAM Fe (TI)
- `systems/EAM-SNAP-SIA-Fe` : <br>
    vacancy with Marinica04 EAM (base) and Marinica MILADY SNAP (TI)

## PAFI Implementations 
! Both examples have *very short samples for testing* : increase `SampleSteps` and `ThermSteps` to 500-1000 + !! !
Run in serial with e.g.
```bash
cd examples/standard
../../build/pafi
```

- `standard` : normal PAFI routine, applied to the Fe SIA case 
- `custom` : example TI application

