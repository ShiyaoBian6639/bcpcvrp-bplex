# instances/cvrplib

CVRPLIB sets A, B, E, F, P (Augerat; Christofides-Eilon; Fisher) downloaded
2026-10-26 from https://galgos.inf.puc-rio.br/cvrplib/ (instance-set ids
4, 5, 9, 10, 12). Each `.vrp` has a `.sol` with the proven optimal routes and
`Cost`. Convention: nint-rounded Euclidean distances; E-n13-k4 and E-n31-k7
carry an EXPLICIT LOWER_ROW matrix.

Other sets used by cvrpbcp already live in the AIRouting tree and are not
copied:

* X (with BKS): `../../../xltw/benchmarks/CVRP`, `../../../HGS-CVRP/Instances/CVRP`
* XML100 with optimal solutions: `../../../data/xml100/{instances,solutions}`
* XML150/200: `../../../exactlabel/inst`
* Solomon and Gehring-Homberger CVRPTW with BKS (one-decimal truncation):
  `../../../xltw/benchmarks/VRPTW`
