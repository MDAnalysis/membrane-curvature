import re
import sys
from pathlib import Path

def clean_pyproject_toml(input_path: str, output_path: str = None) -> None:
    """
    清理 pyproject.toml 文件：
    - 移除 [tool.yapf] 部分
    - 移除 License classifier（PEP 639 已弃用）
    - 用 SPDX ID 替换旧的 PEP 621 license 字段
    - 移除冗余的 Python 3 classifier（只保留 3.11-3.14）
    - 在 [project.urls] 中添加 MC 主页、仓库和 issues 链接

    Args:
        input_path: 输入文件路径
        output_path: 输出文件路径，若为 None 则覆盖原文件
    """
    try:
        with open(input_path, 'r', encoding='utf-8') as f:
            content = f.read()
    except FileNotFoundError:
        print(f"错误：文件 {input_path} 未找到")
        sys.exit(1)
    except Exception as e:
        print(f"读取文件时出错: {e}")
        sys.exit(1)

    # 1. 移除 [tool.yapf] 部分及其内容
    content = re.sub(r'\[tool\.yapf\].*?(?=\[|\Z)', '', content, flags=re.DOTALL)

    # 2. 移除 License classifier（以 "License :: " 开头）
    content = re.sub(r'^\s*"License :: .*",?\s*$', '', content, flags=re.MULTILINE)

    # 3. 替换旧的 PEP 621 license 字段为 SPDX ID
    # 假设旧的格式为 license = {text = "MIT"} 或类似，替换为 license = "MIT"
    content = re.sub(
        r'license\s*=\s*\{[^}]*"([^"]+)"[^}]*\}',
        r'license = "\1"',
        content
    )

    # 4. 移除冗余的 Python 3 classifier，只保留 3.11-3.14
    # 匹配 "Programming Language :: Python :: 3" 以及 3.0-3.10
    def keep_python_classifier(match):
        version = match.group(1)
        # 保留 3.11, 3.12, 3.13, 3.14
        if version in ['3.11', '3.12', '3.13', '3.14']:
            return match.group(0)
        else:
            return ''

    content = re.sub(
        r'"Programming Language :: Python :: (3(?:\.\d+)?)"',
        keep_python_classifier,
        content
    )

    # 5. 在 [project.urls] 中添加链接（如果不存在则创建）
    # 检查是否已有 [project.urls]
    if '[project.urls]' in content:
        # 在 [project.urls] 部分添加缺失的链接
        # 简单处理：在 [project.urls] 下一行添加
        lines = content.split('\n')
        new_lines = []
        in_urls = False
        added_homepage = False
        added_repository = False
        added_issues = False
        for i, line in enumerate(lines):
            new_lines.append(line)
            if line.strip() == '[project.urls]':
                in_urls = True
                continue
            if in_urls:
                # 检查是否已有这些链接
                if 'Homepage' in line:
                    added_homepage = True
                if 'Repository' in line:
                    added_repository = True
                if 'Issues' in line:
                    added_issues = True
                # 遇到下一个 section 或文件结束
                if line.startswith('[') and line != '[project.urls]':
                    # 在进入下一个 section 之前添加缺失的链接
                    if not added_homepage:
                        new_lines.append('Homepage = "https://example.com/mc-homepage"')
                    if not added_repository:
                        new_lines.append('Repository = "https://github.com/example/mc-repo"')
                    if not added_issues:
                        new_lines.append('Issues = "https://github.com/example/mc-repo/issues"')
                    in_urls = False
        # 如果文件结束仍在 [project.urls] 中
        if in_urls:
            if not added_homepage:
                new_lines.append('Homepage = "https://example.com/mc-homepage"')
            if not added_repository:
                new_lines.append('Repository = "https://github.com/example/mc-repo"')
            if not added_issues:
                new_lines.append('Issues = "https://github.com/example/mc-repo/issues"')
        content = '\n'.join(new_lines)
    else:
        # 在文件末尾添加 [project.urls] 部分
        content += '\n[project.urls]\n'
        content += 'Homepage = "https://example.com/mc-homepage"\n'
        content += 'Repository = "https://github.com/example/mc-repo"\n'
        content += 'Issues = "https://github.com/example/mc-repo/issues"\n'

    # 清理多余的空行（可选）
    content = re.sub(r'\n{3,}', '\n\n', content)

    # 写入文件
    output_path = output_path or input_path
    try:
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(content)
        print(f"成功清理并写入 {output_path}")
    except Exception as e:
        print(f"写入文件时出错: {e}")
        sys.exit(1)

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("用法: python cleanup_pyproject.py <input_path> [output_path]")
        sys.exit(1)
    input_path = sys.argv[1]
    output_path = sys.argv[2] if len(sys.argv) > 2 else None
    clean_pyproject_toml(input_path, output_path)
