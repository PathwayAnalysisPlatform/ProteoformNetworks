import pytest


@pytest.fixture(scope="session", autouse=True)
def setupSession():
    print("\nSetup Session")

#
# @pytest.fixture(scope="module", autouse=True)
# def setupModule():
#     print("\nSetup Module")
#
#

@pytest.fixture(scope="function", autouse=True)
def setupFunction():
    print("\nSetup Function")


@pytest.fixture()
def setup1():
    print("\nSetup 1")
    yield
    print("\nTeardown 1")


@pytest.fixture()
def setup2(request):
    print("\nSetup 2")

    def teardown_a():
        print("\nTeardown A")

    def teardown_b():
        print("\nTeardown B")

    request.addfinalizer(teardown_a)
    request.addfinalizer(teardown_b)


def test1(setup1):
    print("Executing test1!")
    assert True


@pytest.mark.usefixtures("setup2")
def test2():
    print("Executing test2!")
    assert True


def test3():
    print("Executing test3!")
    assert True


class TestClass:
    @pytest.fixture(scope="class", autouse=True)
    def setupClass(self):
        print("\nSetup Class")

    @pytest.fixture(scope="function", autouse=False)
    def setupFunction(self):
        print("\nSetup Function in Class")

    def test_it(self):
        print("TestIt")
        assert True

    @pytest.mark.usefixtures("setup2")
    def test_it2(self):
        print("TestIt2")
        assert True
