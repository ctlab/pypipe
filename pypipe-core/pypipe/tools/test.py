from pypipe.tools.toolsconfig import Test
from pypipe.utils import tool


@tool
def one_one():
    return Test.one_one()


@tool
def one_two():
    return Test.one_two()


@tool
def three_one():
    return Test.three_one()

