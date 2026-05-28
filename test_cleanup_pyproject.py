import os
import tempfile
import unittest
from cleanup_pyproject import clean_pyproject_toml

class TestCleanPyproject(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.input_path = os.path.join(self.temp_dir, 'pyproject.toml')
        self.output_path = os.path.join(self.temp_dir, 'pyproject_clean.toml')

    def write_input(self, content):
        with open(self.input_path, 'w', encoding='utf-8') as f:
            f.write(content)

    def read_output(self):
        with open(self.output_path, 'r', encoding='utf-8') as f:
            return f.read()

    def test_remove_yapf(self):
        content = """[tool.yapf]
based_on_style = pep8
indent_width = 4

[project]
name = "test"
"""
        self.write_input(content)
        clean_pyproject_toml(self.input_path, self.output_path)
        output = self.read_output()
        self.assertNotIn('[tool.yapf]', output)

    def test_remove_license_classifier(self):
        content = """[project]
classifiers = [
    "License :: OSI Approved :: MIT License",
    "Programming Language :: Python :: 3",
]
"""
        self.write_input(content)
        clean_pyproject_toml(self.input_path, self.output_path)
        output = self.read_output()
        self.assertNotIn('License ::', output)

    def test_replace_license_field(self):
        content = """[project]
license = "MIT"
"""
        self.write_input(content)
        clean_pyproject_toml(self.input_path, self.output_path)
        output = self.read_output()
        self.assertIn('license = "MIT"', output)

    def test_keep_python_311_314(self):
        content = """[project]
classifiers = [
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.9",
    "Programming Language :: Python :: 3.11",
    "Programming Language :: Python :: 3.12",
    "Programming Language :: Python :: 3.13",
    "Programming Language :: Python :: 3.14",
]
"""
        self.write_input(content)
        clean_pyproject_toml(self.input_path, self.output_path)
        output = self.read_output()
        self.assertIn('Python :: 3.11', output)
        self.assertIn('Python :: 3.12', output)
        self.assertIn('Python :: 3.13', output)
        self.assertIn('Python :: 3.14', output)
        self.assertNotIn('Python :: 3.9', output)
        self.assertNotIn('Python :: 3"', output)  # 确保移除了 "Python :: 3"

    def test_add_urls(self):
        content = """[project]
name = "test"
"""
        self.write_input(content)
        clean_pyproject_toml(self.input_path, self.output_path)
        output = self.read_output()
        self.assertIn('[project.urls]', output)
        self.assertIn('Homepage', output)
        self.assertIn('Repository', output)
        self.assertIn('Issues', output)

if __name__ == '__main__':
    unittest.main()
