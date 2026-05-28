import tempfile
import os
import sys
sys.path.insert(0, os.path.dirname(__file__))
from cleanup_pyproject import clean_pyproject_toml

def test_clean_pyproject_toml():
    # 创建临时输入文件
    input_content = """[build-system]
requires = ["setuptools"]
build-backend = "setuptools.build_meta"

[project]
name = "example"
version = "0.1.0"
license = {text = "MIT"}
classifiers = [
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3 :: Only",
    "Programming Language :: Python :: 3.10",
    "Programming Language :: Python :: 3.11",
    "Programming Language :: Python :: 3.12",
    "License :: OSI Approved :: MIT License",
]

[tool.yapf]
based_on_style = "pep8"
"""
    expected_output = """[build-system]
requires = ["setuptools"]
build-backend = "setuptools.build_meta"

[project]
name = "example"
version = "0.1.0"
license = "MIT"
classifiers = [
    "Programming Language :: Python :: 3.11",
    "Programming Language :: Python :: 3.12",
]

[project.urls]
Homepage = "https://example.com"
Repository = "https://github.com/example/repo"
Issues = "https://github.com/example/repo/issues"
"""
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.toml') as f:
        f.write(input_content)
        input_path = f.name
    output_path = input_path + '.out'
    try:
        clean_pyproject_toml(input_path, output_path)
        with open(output_path, 'r') as f:
            result = f.read()
        # 比较时忽略末尾换行差异
        assert result.strip() == expected_output.strip(), f"结果不匹配\n期望：\n{expected_output}\n实际：\n{result}"
        print("测试通过！")
    finally:
        os.unlink(input_path)
        if os.path.exists(output_path):
            os.unlink(output_path)

if __name__ == '__main__':
    test_clean_pyproject_toml()
