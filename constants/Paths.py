from pathlib import Path


# todo rname path or Dart or somethin
class project_path:
    def __init__(self):
        self.cwd = Path.cwd()
        self.project_path = None

        for parent_path in self.cwd.parents:
            files_on_same_level = list(parent_path.glob('*'))
            for file_on_same_level in files_on_same_level:
                if file_on_same_level.is_dir():
                    if str(file_on_same_level.name) == ".git":
                        self.project_path = file_on_same_level.parent
                    else:
                        pass
                else:
                    pass
        assert self.project_path is not None

    def extend(self, *relpaths: object) -> object:
        return Path(self.project_path, *relpaths)


project_path()
#print(project_path().extend("data","Filtered_Jsons","INTEGRATION_TEST_LIGAND_DATABASE_270623.json"))
