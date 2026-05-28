import tempfile
import os
import sys
from pathlib import Path

# 将主脚本所在目录加入路径
sys.path.insert(0, os.path.dirname(__file__))

from cleanup_pyproject import clean_pyproject_toml

def test_clean_pyproject_toml():
    # 创建临时目录
    with tempfile.TemporaryDirectory() as tmpdir:
        input_path = os.path.join(tmpdir, 'pyproject.toml')
        output_path = os.path.join(tmpdir, 'pyproject_clean.toml')

        # 编写测试用的 pyproject.toml 内容
        test_content = """[build-system]
requires = ["setuptools"]
build-backend = "setuptools.build_meta"

[project]
name = "example"
version = "0.1.0"
license = {text = "MIT"}
classifiers = [
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.9",
    "Programming Language :: Python :: 3.10",
    "Programming Language :: Python :: 3.11",
    "Programming Language :: Python :: 3.12",
    "License :: OSI Approved :: MIT License",
]

[tool.yapf]
based_on_style = "pep8"

[project.urls]
Documentation = "https://example.com/docs"
"""

        with open(input_path, 'w', encoding='utf-8') as f:
            f.write(test_content)

        # 执行清理
        clean_pyproject_toml(input_path, output_path)

        # 读取清理后的内容
        with open(output_path, 'r', encoding='utf-8') as f:
            cleaned = f.read()

        # 验证
        assert '[tool.yapf]' not in cleaned, "应移除 [tool.yapf]"
        assert 'License :: OSI Approved :: MIT License' not in cleaned, "应移除 License classifier"
        assert 'license = "MIT"' in cleaned, "license 应替换为 SPDX ID"
        assert 'Programming Language :: Python :: 3.9' not in cleaned, "应移除 3.9 classifier"
        assert 'Programming Language :: Python :: 3.11' in cleaned, "应保留 3.11 classifier"
        assert 'Programming Language :: Python :: 3.12' in cleaned, "应保留 3.12 classifier"
        assert 'Homepage = "https://example.com/mc-homepage"' in cleaned, "应添加 Homepage"
        assert 'Repository = "https://github.com/example/mc-repo"' in cleaned, "应添加 Repository"
        assert 'Issues = "https://github.com/example/mc-repo/issues"' in cleaned, "应添加 Issues"
        print("所有测试通过！")

if __name__ == '__main__':
    test_clean_pyproject_toml()
