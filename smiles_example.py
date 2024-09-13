from rdkit import Chem
m = Chem.MolFromSmiles('c1nccc1O(C)C',sanitize=False)
errs = Chem.DetectChemistryProblems(m)
for er in errs:
    if er.GetType()=='AtomValenceException':
        print(er.GetType(),er.Message(),er.GetAtomIdx())
    elif er.GetType()=='KekulizeException':
        print(er.GetType(),er.Message(),er.GetAtomIndices())
