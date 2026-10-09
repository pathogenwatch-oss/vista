from vista.vista import search_parallelism


def test_search_parallelism_stays_within_cpu_budget() -> None:
    assert search_parallelism(cpus=16, library_count=3) == (3, 5)
    assert search_parallelism(cpus=2, library_count=3) == (2, 1)
    assert search_parallelism(cpus=4, library_count=1) == (1, 4)
