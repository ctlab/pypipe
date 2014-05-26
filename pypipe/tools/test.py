from pypipe.utils import tool
from pypipe.tools.toolsconfig import Test


@tool
def one_one():
    return Test.one_one()


@tool
def one_two():
    return Test.one_two()


@tool
def three_one():
    return Test.three_one()

