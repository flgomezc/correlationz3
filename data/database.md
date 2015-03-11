DATABASE
========

CosmoSim is a new database holding data from Bolshoi, MDR1, MDPL and recently
added BolshoiP with Planck cosmology.

The Snapshot 87 is equivalent to Redshift z~3 

```SQL
SELECT bdmID, hostFlag, x,y,z, Mvir
FROM BolshoiP.BDMV
WHERE snapnum=87
```
