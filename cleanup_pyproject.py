import re
import sys

def clean_pyproject_toml(input_path: str, output_path: str = None) -> None:
    """
    清理 pyproject.toml 文件：
    - 移除 [tool.yapf] 部分
    - 移除 License classifier（PEP 639 已弃用）
    - 用 SPDX ID 替换已弃用的 license 字段
    - 移除冗余的 Python 3 classifier（保留 3.11-3.14）
    - 在 [project.urls] 中添加 homepage、repository、issues 链接
    """
    try:
        with open(input_path, 'r', encoding='utf-8') as f:
            content = f.read()
    except FileNotFoundError:
        print(f"错误：文件 {input_path} 未找到")
        sys.exit(1)
    except Exception as e:
        print(f"读取文件时出错：{e}")
        sys.exit(1)

    # 1. 移除 [tool.yapf] 部分及其内容
    content = re.sub(r'\[tool\.yapf\].*?(?=\n\[|\Z)', '', content, flags=re.DOTALL)

    # 2. 移除 License classifier（例如："License :: OSI Approved :: MIT License"）
    content = re.sub(r'\s*"License :: [^"]+"\s*,?\n?', '', content)

    # 3. 替换已弃用的 license 字段为 SPDX ID（假设原始 license 字段包含 "MIT"）
    #    匹配 license = {text = "..."} 或 license = "..." 并替换
    content = re.sub(
        r'license\s*=\s*\{[^}]*text\s*=\s*"([^"]+)"[^}]*\}',
        r'license = "\1"',
        content
    )
    # 如果 license 字段已经是字符串，确保它是 SPDX 格式（假设原始值正确）
    # 这里不做进一步验证，仅确保格式

    # 4. 移除冗余的 Python 3 classifier（保留 3.11-3.14）
    #    删除所有 "Programming Language :: Python :: 3" 和 "Programming Language :: Python :: 3 :: Only"
    #    以及 3.10 及以下的版本
    classifiers_to_remove = [
        r'"Programming Language :: Python :: 3"',
        r'"Programming Language :: Python :: 3 :: Only"',
        r'"Programming Language :: Python :: 3\.(?:0|1|2|3|4|5|6|7|8|9|10)"',
    ]
    for pattern in classifiers_to_remove:
        content = re.sub(pattern + r'\s*,?\n?', '', content)

    # 5. 在 [project.urls] 中添加链接
    #    如果 [project.urls] 不存在，则创建；如果存在，则追加
    urls_section = "[project.urls]\nHomepage = \"https://example.com\"\nRepository = \"https://github.com/example/repo\"\nIssues = \"https://github.com/example/repo/issues\"\n"
    if '[project.urls]' not in content:
        content += '\n' + urls_section
    else:
        # 在 [project.urls] 部分末尾添加，避免重复
        # 简单处理：如果已有 Homepage 等，不重复添加
        if 'Homepage' not in content:
            content = re.sub(r'(\[project\.urls\]\n)', r'\1Homepage = "https://example.com"\n', content)
        if 'Repository' not in content:
            content = re.sub(r'(\[project\.urls\]\n)', r'\1Repository = "https://github.com/example/repo"\n', content)
        if 'Issues' not in content:
            content = re.sub(r'(\[project\.urls\]\n)', r'\1Issues = "https://github.com/example/repo/issues"\n', content)

    # 清理多余的空行
    content = re.sub(r'\n{3,}', '\n\n', content)
    content = content.strip() + '\n'

    # 写入输出文件
    if output_path is None:
        output_path = input_path
    try:
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(content)
        print(f"成功：已清理并保存到 {output_path}")
    except Exception as e:
        print(f"写入文件时出错：{e}")
        sys.exit(1)

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("用法：python cleanup_pyproject.py <输入文件> [输出文件]")
        sys.exit(1)
    input_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else None
    clean_pyproject_toml(input_file, output_file)
