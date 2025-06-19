# 🧩 DART Assembly Refactor — Progress & Notes

This document tracks the ongoing refactor of DART's assembly logic, aiming for clearer structure, modularity, and better long-term maintainability.

---

## ✅ TODO Progress Tracker

| Task Description                                                                                                                                        | Status      | Priority |
|---------------------------------------------------------------------------------------------------------------------------------------------------------|-------------|----------|
| Implement warnings that are carried to the end                                                                                                          | To Do       | High     |
| Document assembly process                                                                                                                               | To Do       | Medium   |
| Refactor IsomerFactory.generate() function as the consistnency between nested sets of ligands for each isomer is confusing                              | To Do       | Low      |
| Change DARTIsomer to AssembledIsomer                                                                                                                    | In Progress | High     |
| Update arrays in atoms objects to be more descriptive                                                                                                   | To Do       | High     |
| Refactor BaseMolecule so that I can use its atoms object rather than having to define self.DART_atoms (could introduce specific AssembledIsomer labels) | To Do       | High     |
| Refactor AssemblySave to be able to deal with the class method                                                                                          | To Do       | High     |

---

## 🧠 Notes & Overview

- **Modifications of code external to DART assembly workflow**: 
  - I modified assert_graph_and_coordinates_are_consistent in utils_graph.py such that the comparison
  in the second assert statement is done in an ordered way


---