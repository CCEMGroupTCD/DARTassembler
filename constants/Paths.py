from pathlib import Path
#todo rname path or Dart or somethin
class project_path:
    def __init__(self):
        self.cwd = Path.cwd()
        self.project_path = self.cwd.parent

    def extend(self, *relpaths: object) -> object:
        return Path(self.project_path, *relpaths)

print(project_path().extend("data","Filtered_Jsons","INTEGRATION_TEST_LIGAND_DATABASE_270623.json"))