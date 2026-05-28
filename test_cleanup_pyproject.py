import unittest
import tempfile
import os
import toml
from cleanup_pyproject import cleanup_pyproject_toml

class TestCleanupPyproject(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.input_path = os.path.join(self.temp_dir, 'pyproject.toml')
        self.output_path = os.path.join(self.temp_dir, 'output.toml')

    def tearDown(self):
        import shutil
        shutil.rmtree(self.temp_dir)

    def _write_input(self, content: str):
        with open(self.input_path, 'w', encoding='utf-8') as f:
            f.write(content)

    def test_remove_yapf(self):
        content = '''
[project]
name = "test"

[tool.yapf]
based_on_style = "pep8"
'''
        self._write_input(content)
        result = cleanup_pyproject_toml(self.input_path, self.output_path)
        self.assertNotIn('yapf', result.get('tool', {}))

    def test_remove_license_classifier(self):
        content = '''
[project]
classifiers = [
    "License :: OSI Approved :: MIT License",
    "Programming Language :: Python :: 3.11",
]
'''
        self._write_input(content)
        result = cleanup_pyproject_toml(self.input_path, self.output_path)
        classifiers = result['project'].get('classifiers', [])
        self.assertNotIn("License :: OSI Approved :: MIT License", classifiers)
        self.assertIn("Programming Language :: Python :: 3.11", classifiers)

    def test_replace_license_with_spdx(self):
        content = '''
[project]
license = {text = "MIT License"}
'''
        self._write_input(content)
        result = cleanup_pyproject_toml(self.input_path, self.output_path)
        self.assertEqual(result['project']['license'], 'MIT')

    def test_keep_python_classifiers_311_314(self):
        content = '''
[project]
classifiers = [
    "Programming Language :: Python :: 3.9",
    "Programming Language :: Python :: 3.10",
    "Programming Language :: Python :: 3.11",
    "Programming Language :: Python :: 3.12",
    "Programming Language :: Python :: 3.13",
    "Programming Language :: Python :: 3.14",
]
'''
        self._write_input(content)
        result = cleanup_pyproject_toml(self.input_path, self.output_path)
        classifiers = result['project']['classifiers']
        self.assertNotIn("Programming Language :: Python :: 3.9", classifiers)
        self.assertNotIn("Programming Language :: Python :: 3.10", classifiers)
        self.assertIn("Programming Language :: Python :: 3.11", classifiers)
        self.assertIn("Programming Language :: Python :: 3.14", classifiers)

    def test_add_urls(self):
        content = '''
[project]
name = "test"
'''
        self._write_input(content)
        result = cleanup_pyproject_toml(self.input_path, self.output_path)
        urls = result['project']['urls']
        self.assertIn('homepage', urls)
        self.assertIn('repository', urls)
        self.assertIn('issues', urls)

    def test_no_overwrite_existing_urls(self):
        content = '''
[project]
urls = {homepage = "https://existing.com"}
'''
        self._write_input(content)
        result = cleanup_pyproject_toml(self.input_path, self.output_path)
        self.assertEqual(result['project']['urls']['homepage'], 'https://existing.com')

if __name__ == '__main__':
    unittest.main()