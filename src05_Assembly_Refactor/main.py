from src05_Assembly_Refactor.DART_Assembly import DARTAssembly
from constants.Paths import project_path
from src05_Assembly_Refactor.utilities_assembly import visualize

def main(assembly_input_path: str):
    assembly = DARTAssembly(assembly_input_path=USER_INPUT)
    assembly.run_all_batches()

    return assembly

if __name__ == "__main__":

    USER_INPUT = 'assembly_input.yml'   # In this file the user specifies which input they want
    assembly = main(assembly_input_path=USER_INPUT)




    #%% ==============    Doublecheck refactoring    ==================
    print('\nDoublechecking if output is same as before:')
    from src14_Assembly_Unit_Test.Assembly_test import AssemblyIntegrationTest

    # Check that the output movie is the same
    if assembly.optimization_movie:
        print('Checking if movie is the same...')
        benchmark_movie_path = project_path().extend("tmp", "opt_movie_Benchmark.xyz")
        test_movie_path = assembly.optimization_movie_path
        df_movie_diff = AssemblyIntegrationTest(benchmark_movie_path, test_movie_path).compare_xyz_files()

    print('Checking if output xyz files are the same...')
    benchmark_file = '/Users/timosommer/PhD/projects/RCA/projects/DART/src14_Assembly_Unit_Test/INTEGRATION_TEST_after_refactoring_input.xyz'
    test_file = assembly.concat_xyz_path
    allowed_differences = 1e-5
    df_xzy_diff = AssemblyIntegrationTest(benchmark_file, test_file, tol=allowed_differences).compare_xyz_files()

    print('Done!')
