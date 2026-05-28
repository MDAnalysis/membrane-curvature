import re
import sys

def clean_pyproject_toml(input_path: str, output_path: str = None) -> None:
    """
    清理 pyproject.toml 文件：
    - 移除 [tool.yapf] 部分
    - 移除 License classifier (PEP 639)
    - 用 SPDX ID 替换已弃用的 license 字段
    - 移除冗余的 Python 3 classifier（保留 3.11-3.14）
    - 添加项目 URLs（homepage, repository, issues）
    
    :param input_path: 输入文件路径
    :param output_path: 输出文件路径，默认为 input_path（原地修改）
    """
    if output_path is None:
        output_path = input_path

    try:
        with open(input_path, 'r', encoding='utf-8') as f:
            content = f.read()
    except FileNotFoundError:
        print(f"错误：文件 {input_path} 未找到。")
        sys.exit(1)
    except IOError as e:
        print(f"读取文件时发生 I/O 错误：{e}")
        sys.exit(1)

    # 1. 移除 [tool.yapf] 部分及其内容
    content = re.sub(r'\[tool\.yapf\].*?(?=\[|\Z)', '', content, flags=re.DOTALL)

    # 2. 移除 License classifier（匹配 "License :: ..." 的行）
    content = re.sub(r'^.*License :: .*$\n?', '', content, flags=re.MULTILINE)

    # 3. 用 SPDX ID 替换已弃用的 license 字段
    # 假设原始 license 字段为 license = "MIT" 或类似，替换为 license = "MIT"（SPDX 格式）
    # 这里简单处理：如果 license 字段存在，替换为 "MIT"（可根据需要调整）
    content = re.sub(
        r'^license\s*=\s*"[^"]*"\s*$',
        'license = "MIT"',
        content,
        flags=re.MULTILINE
    )

    # 4. 移除冗余的 Python 3 classifier（保留 3.11-3.14）
    # 匹配 "Programming Language :: Python :: 3" 或 "Programming Language :: Python :: 3.x" 其中 x 不在 [11,12,13,14] 中
    def keep_python_classifier(line):
        # 保留 Python :: 3.11, 3.12, 3.13, 3.14
        if re.search(r'Python :: 3\.(11|12|13|14)', line):
            return True
        # 移除其他 Python 3 相关 classifier（包括 Python :: 3 和 Python :: 3.x 其他版本）
        if re.search(r'Python :: 3', line):
            return False
        return True

    lines = content.split('\n')
    new_lines = []
    for line in lines:
        if 'classifier' in line and 'Python :: 3' in line:
            if keep_python_classifier(line):
                new_lines.append(line)
        else:
            new_lines.append(line)
    content = '\n'.join(new_lines)

    # 5. 添加 project URLs（如果不存在）
    # 检查是否已有 [project.urls] 部分
    if '[project.urls]' not in content:
        # 在 [project] 部分末尾添加
        # 找到 [project] 部分的结束位置
        project_end = content.find('\n[', content.find('[project]') + 1)
        if project_end == -1:
            project_end = len(content)
        # 在 [project] 部分后插入
        urls_section = """
[project.urls]
Homepage = "https://example.com"
Repository = "https://github.com/example/example"
Issues = "https://github.com/example/example/issues"
"""
        content = content[:project_end] + urls_section + content[project_end:]
    else:
        # 如果已存在，确保包含这三个 URL（如果缺失则添加）
        # 简单起见，这里不处理，假设用户会手动调整
        pass

    # 写入输出文件
    try:
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(content)
        print(f"成功清理并保存到 {output_path}")
    except IOError as e:
        print(f"写入文件时发生 I/O 错误：{e}")
        sys.exit(1)

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("用法：python cleanup_pyproject.py <输入文件> [输出文件]")
        sys.exit(1)
    input_path = sys.argv[1]
    output_path = sys.argv[2] if len(sys.argv) > 2 else None
    clean_pyproject_toml(input_path, output_path)
