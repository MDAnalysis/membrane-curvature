import re
import sys

def cleanup_pyproject_toml(input_path: str, output_path: str = None) -> None:
    """
    清理 pyproject.toml 文件：
    - 移除 [tool.yapf] 部分
    - 移除 License classifier（PEP 639 废弃）
    - 替换废弃的 PEP 621 license 字段为 SPDX ID
    - 移除冗余的 Python 3 classifier（保留 3.11-3.14）
    - 在 [project.urls] 中添加 MC 主页、仓库和问题链接

    Args:
        input_path: 输入文件路径
        output_path: 输出文件路径，默认为覆盖原文件
    """
    if output_path is None:
        output_path = input_path

    try:
        with open(input_path, 'r', encoding='utf-8') as f:
            content = f.read()
    except FileNotFoundError:
        print(f"错误：文件 {input_path} 未找到")
        sys.exit(1)
    except IOError as e:
        print(f"读取文件错误：{e}")
        sys.exit(1)

    # 1. 移除 [tool.yapf] 部分及其内容
    content = re.sub(r'\[tool\.yapf\].*?(?=\[|\Z)', '', content, flags=re.DOTALL)

    # 2. 移除 License classifier（例如："License :: OSI Approved :: MIT License"）
    content = re.sub(r'\s*"License :: .*?",?\n?', '', content)

    # 3. 替换废弃的 PEP 621 license 字段为 SPDX ID
    # 假设原字段为 license = {text = "MIT"} 或 license = "MIT"，替换为 license = "MIT"
    content = re.sub(
        r'license\s*=\s*\{[^}]*text\s*=\s*"([^"]+)"[^}]*\}',
        r'license = "\1"',
        content
    )

    # 4. 移除冗余的 Python 3 classifier（保留 3.11-3.14）
    # 保留包含 "3.11", "3.12", "3.13", "3.14" 的 classifier，移除其他 Python 3 版本
    def filter_classifier(match):
        line = match.group(0)
        # 如果包含 3.11-3.14 则保留
        if re.search(r'"Programming Language :: Python :: 3\.(11|12|13|14)"', line):
            return line
        # 如果是 "Python :: 3" 或 "Python :: 3.x"（x 不在 11-14）则移除
        if re.search(r'"Programming Language :: Python :: 3"', line) and not re.search(r'3\.(11|12|13|14)', line):
            return ''
        return line

    # 逐行处理 classifiers 部分
    lines = content.split('\n')
    new_lines = []
    in_classifiers = False
    for line in lines:
        if 'classifiers = [' in line:
            in_classifiers = True
            new_lines.append(line)
        elif in_classifiers:
            if ']' in line:
                in_classifiers = False
                new_lines.append(line)
            else:
                filtered = filter_classifier(line)
                if filtered:
                    new_lines.append(filtered)
        else:
            new_lines.append(line)
    content = '\n'.join(new_lines)

    # 5. 在 [project.urls] 中添加链接
    # 如果存在 [project.urls] 部分，在末尾添加；否则创建
    urls_section = """[project.urls]
Homepage = "https://example.com/mc"
Repository = "https://github.com/example/mc"
Issues = "https://github.com/example/mc/issues"
"""
    if '[project.urls]' in content:
        # 替换已有的 [project.urls] 部分
        content = re.sub(
            r'\[project\.urls\].*?(?=\[|\Z)',
            urls_section.rstrip(),
            content,
            flags=re.DOTALL
        )
    else:
        # 在文件末尾添加
        content += '\n' + urls_section

    # 写入输出文件
    try:
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(content)
        print(f"成功清理并写入 {output_path}")
    except IOError as e:
        print(f"写入文件错误：{e}")
        sys.exit(1)

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("用法：python cleanup_pyproject.py <输入文件> [输出文件]")
        sys.exit(1)
    input_path = sys.argv[1]
    output_path = sys.argv[2] if len(sys.argv) > 2 else None
    cleanup_pyproject_toml(input_path, output_path)
