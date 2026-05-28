import tempfile
import os
import sys
sys.path.insert(0, os.path.dirname(__file__))
from cleanup_pyproject import cleanup_pyproject_toml

def test_cleanup():
    # 创建临时测试文件
    test_input = """[build-system]
requires = ["setuptools"]
build-backend = "setuptools.build_meta"

[project]
name = "mc"
version = "0.1.0"
license = {text = "MIT"}
classifiers = [
    "License :: OSI Approved :: MIT License",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.9",
    "Programming Language :: Python :: 3.10",
    "Programming Language :: Python :: 3.11",
    "Programming Language :: Python :: 3.12",
    "Programming Language :: Python :: 3.13",
    "Programming Language :: Python :: 3.14",
]

[tool.yapf]
based_on_style = "pep8"
"""

    expected_output = """[build-system]
requires = ["setuptools"]
build-backend = "setuptools.build_meta"

[project]
name = "mc"
version = "0.1.0"
license = "MIT"
classifiers = [
    "Programming Language :: Python :: 3.11",
    "Programming Language :: Python :: 3.12",
    "Programming Language :: Python :: 3.13",
    "Programming Language :: Python :: 3.14",
]

[project.urls]
Homepage = "https://example.com/mc"
Repository = "https://github.com/example/mc"
Issues = "https://github.com/example/mc/issues"
"""

    with tempfile.NamedTemporaryFile(mode='w', suffix='.toml', delete=False) as f:
        f.write(test_input)
        input_path = f.name

    output_path = input_path + '.out'
    cleanup_pyproject_toml(input_path, output_path)

    with open(output_path, 'r') as f:
        result = f.read()

    # 清理多余空行以便比较
    result_clean = '\n'.join(line for line in result.split('\n') if line.strip())
    expected_clean = '\n'.join(line for line in expected_output.split('\n') if line.strip())

    assert result_clean == expected_clean, f"结果不匹配\n结果:\n{result_clean}\n期望:\n{expected_clean}"
    print("测试通过！")

    # 清理临时文件
    os.unlink(input_path)
    os.unlink(output_path)

if __name__ == '__main__':
    test_cleanup()
