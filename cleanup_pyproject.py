import toml
import re

def cleanup_pyproject_toml(input_path: str, output_path: str = None) -> dict:
    """
    清理 pyproject.toml 文件：
    - 移除 [tool.yapf] 部分
    - 移除 License classifier（PEP 639 已弃用）
    - 将旧的 PEP 621 license 字段替换为 SPDX ID
    - 移除冗余的 Python 3 分类器（保留 3.11-3.14）
    - 添加项目 URL（homepage, repository, issues）
    
    :param input_path: 输入文件路径
    :param output_path: 输出文件路径，若为 None 则覆盖原文件
    :return: 清理后的 pyproject.toml 字典
    """
    try:
        with open(input_path, 'r', encoding='utf-8') as f:
            content = f.read()
        data = toml.loads(content)
    except FileNotFoundError:
        raise FileNotFoundError(f"文件 {input_path} 未找到")
    except toml.TomlDecodeError as e:
        raise ValueError(f"TOML 解析错误: {e}")

    # 1. 移除 [tool.yapf]
    if 'tool' in data and 'yapf' in data['tool']:
        del data['tool']['yapf']
        if not data['tool']:
            del data['tool']

    # 2. 移除 License classifier
    if 'project' in data and 'classifiers' in data['project']:
        classifiers = data['project']['classifiers']
        # 匹配 License :: 开头的分类器
        data['project']['classifiers'] = [
            c for c in classifiers if not c.startswith('License ::')
        ]
        if not data['project']['classifiers']:
            del data['project']['classifiers']

    # 3. 替换旧的 PEP 621 license 字段为 SPDX ID
    if 'project' in data and 'license' in data['project']:
        license_field = data['project']['license']
        # 如果 license 是字典且包含 text 字段，尝试提取 SPDX ID
        if isinstance(license_field, dict):
            text = license_field.get('text', '')
            # 常见 SPDX ID 映射（可扩展）
            spdx_map = {
                'MIT': 'MIT',
                'Apache-2.0': 'Apache-2.0',
                'GPL-3.0-only': 'GPL-3.0-only',
                'BSD-3-Clause': 'BSD-3-Clause',
            }
            # 尝试从 text 中匹配已知的 SPDX ID
            matched = None
            for key, spdx in spdx_map.items():
                if key.lower() in text.lower():
                    matched = spdx
                    break
            if matched:
                data['project']['license'] = matched
            else:
                # 如果无法匹配，保留原字段但发出警告
                print(f"警告: 无法从 license 字段 '{text}' 提取 SPDX ID，保留原值")
        # 如果已经是字符串，则视为 SPDX ID，不做修改

    # 4. 移除冗余的 Python 3 分类器（保留 3.11-3.14）
    if 'project' in data and 'classifiers' in data['project']:
        classifiers = data['project']['classifiers']
        # 保留非 Python 版本分类器以及 3.11-3.14
        keep_pattern = re.compile(r'Programming Language :: Python :: 3\.(?:11|12|13|14)$')
        data['project']['classifiers'] = [
            c for c in classifiers
            if not c.startswith('Programming Language :: Python :: 3') or keep_pattern.match(c)
        ]
        if not data['project']['classifiers']:
            del data['project']['classifiers']

    # 5. 添加项目 URL
    if 'project' not in data:
        data['project'] = {}
    if 'urls' not in data['project']:
        data['project']['urls'] = {}
    # 仅当不存在时添加，避免覆盖
    if 'homepage' not in data['project']['urls']:
        data['project']['urls']['homepage'] = 'https://example.com'
    if 'repository' not in data['project']['urls']:
        data['project']['urls']['repository'] = 'https://github.com/example/repo'
    if 'issues' not in data['project']['urls']:
        data['project']['urls']['issues'] = 'https://github.com/example/repo/issues'

    # 写入文件
    output_path = output_path or input_path
    try:
        with open(output_path, 'w', encoding='utf-8') as f:
            toml.dump(data, f)
    except Exception as e:
        raise IOError(f"写入文件 {output_path} 失败: {e}")

    return data

if __name__ == '__main__':
    import sys
    if len(sys.argv) < 2:
        print("用法: python cleanup_pyproject.py <input_path> [output_path]")
        sys.exit(1)
    input_path = sys.argv[1]
    output_path = sys.argv[2] if len(sys.argv) > 2 else None
    try:
        result = cleanup_pyproject_toml(input_path, output_path)
        print("清理完成，输出文件:", output_path or input_path)
    except Exception as e:
        print(f"错误: {e}")
        sys.exit(1)