# Changelog

## [0.7.6](https://github.com/varfish-org/mehari-data-tx/compare/v0.7.5...v0.7.6) (2024-11-07)


### Bug Fixes

* don't use anaconda defaults channel ([#51](https://github.com/varfish-org/mehari-data-tx/issues/51)) ([1cfd71a](https://github.com/varfish-org/mehari-data-tx/commit/1cfd71ac0f73fc13b08aaa8791b0503ffbf75c0d))
* re-enable refseq MT grafting from ensembl ([#54](https://github.com/varfish-org/mehari-data-tx/issues/54)) ([6b6b5d2](https://github.com/varfish-org/mehari-data-tx/commit/6b6b5d2b8bc7ea5d1ac2746b701a3ace6f4fa923))
* update hgnc download url ([#52](https://github.com/varfish-org/mehari-data-tx/issues/52)) ([5e3a89a](https://github.com/varfish-org/mehari-data-tx/commit/5e3a89abb8dbbf690d5defe7150feeef80368739))

## [0.7.5](https://github.com/varfish-org/mehari-data-tx/compare/v0.7.4...v0.7.5) (2024-08-16)


### Bug Fixes

* try without explicit shell and cp instead of mv ([6a963a3](https://github.com/varfish-org/mehari-data-tx/commit/6a963a38bdfa7fa9ac8e3e289d58fd41f84f5080))

## [0.7.4](https://github.com/varfish-org/mehari-data-tx/compare/v0.7.3...v0.7.4) (2024-08-16)


### Miscellaneous Chores

* release v0.7.4 ([ae6899c](https://github.com/varfish-org/mehari-data-tx/commit/ae6899c93e6ac373d32c82fe39a3e31cf54140c1))

## [0.7.3](https://github.com/varfish-org/mehari-data-tx/compare/v0.7.3...v0.7.3) (2024-08-16)


### Features

* bump cdot to v0.2.24 ([#33](https://github.com/varfish-org/mehari-data-tx/issues/33)) ([c08ae15](https://github.com/varfish-org/mehari-data-tx/commit/c08ae15c7580000dfdee493bbcafc4d16ff6cf80))
* bump mehari to 0.2.1 ([#5](https://github.com/varfish-org/mehari-data-tx/issues/5)) ([d3caf79](https://github.com/varfish-org/mehari-data-tx/commit/d3caf7992a4f7d2b169914536cff5b00ff0811f4))
* bump mehari to v0.21.0 - release as v0.4.0 ([#15](https://github.com/varfish-org/mehari-data-tx/issues/15)) ([07e5ead](https://github.com/varfish-org/mehari-data-tx/commit/07e5ead23b5137b210a85b913f58ac62a3383d7b))
* Convert build scripts to snakemake workflow ([#38](https://github.com/varfish-org/mehari-data-tx/issues/38)) ([45c7021](https://github.com/varfish-org/mehari-data-tx/commit/45c702113141e7009bc882185d1e58372c54f4ac))
* implement tx protobuf file building ([#2](https://github.com/varfish-org/mehari-data-tx/issues/2)) ([44d5495](https://github.com/varfish-org/mehari-data-tx/commit/44d549596ef66fc7303529e2b6077fc6b508ba74))
* use "seqrepo fetch-load" to fetch missing sequences ([#32](https://github.com/varfish-org/mehari-data-tx/issues/32)) ([7b8c063](https://github.com/varfish-org/mehari-data-tx/commit/7b8c063763bce26faea9d069466f8094cb6f5731))
* use mehari ghcr.io docker image and bump to v0.25.1 ([#30](https://github.com/varfish-org/mehari-data-tx/issues/30)) ([e0fc1a2](https://github.com/varfish-org/mehari-data-tx/commit/e0fc1a217ec17471e8720bfa325814e01727b9c5))


### Bug Fixes

* adding mitochondrial transcripts from ENSEMBL ([#28](https://github.com/varfish-org/mehari-data-tx/issues/28)) ([5bcf163](https://github.com/varfish-org/mehari-data-tx/commit/5bcf1633226c2e932f52d99e7c8381a817777942))
* bump cdot to 0.2.24 for chrMT transcripts ([#31](https://github.com/varfish-org/mehari-data-tx/issues/31)) ([d40bd32](https://github.com/varfish-org/mehari-data-tx/commit/d40bd32296a761be058953392647217872496412))
* bump cdot to v0.2.25 ([#37](https://github.com/varfish-org/mehari-data-tx/issues/37)) ([9cd243c](https://github.com/varfish-org/mehari-data-tx/commit/9cd243cb87a98dd7280dfc823bedd60d144b0866))
* bump mehari to 0.25.4 ([#35](https://github.com/varfish-org/mehari-data-tx/issues/35)) ([253b6ac](https://github.com/varfish-org/mehari-data-tx/commit/253b6acdf7d3c6541e7b146f9179a73b91d7f5f6))
* correct conditional arguments in build script ([#26](https://github.com/varfish-org/mehari-data-tx/issues/26)) ([31a1159](https://github.com/varfish-org/mehari-data-tx/commit/31a1159c4a42adb39167ca4f69c3f5f3d077b69b))
* mehari version to 0.21.0 on release ([#17](https://github.com/varfish-org/mehari-data-tx/issues/17)) ([c0e45f2](https://github.com/varfish-org/mehari-data-tx/commit/c0e45f24fc7732014c920ff94a84f3b686c78ffd))
* properly use GENOME_RELEASE to fix MANE transcripts ([#24](https://github.com/varfish-org/mehari-data-tx/issues/24)) ([79e1854](https://github.com/varfish-org/mehari-data-tx/commit/79e18542de40057d6eeeefbe9bdff208e1b44df3))
* release asset src path, list files during build_data_release ([#46](https://github.com/varfish-org/mehari-data-tx/issues/46)) ([a0dfb9b](https://github.com/varfish-org/mehari-data-tx/commit/a0dfb9bc4088bfc51a8d8e1b339c8f08292c9a1d))
* release please src path ([#40](https://github.com/varfish-org/mehari-data-tx/issues/40)) ([adcba0f](https://github.com/varfish-org/mehari-data-tx/commit/adcba0fa4ea38222c919ea9197c84d8ca23b5187))
* release please src path, part 2 ([#44](https://github.com/varfish-org/mehari-data-tx/issues/44)) ([6437ef6](https://github.com/varfish-org/mehari-data-tx/commit/6437ef6de74c3735b1cd2c8b1eb96d3bf916cbb4))
* request sha256sums for txdb and report ([#42](https://github.com/varfish-org/mehari-data-tx/issues/42)) ([cfdfb24](https://github.com/varfish-org/mehari-data-tx/commit/cfdfb24926be12409b365d84acfbe61a2c9076bf))
* using mane tx list in phase 2 as well ([#20](https://github.com/varfish-org/mehari-data-tx/issues/20)) ([34eeefc](https://github.com/varfish-org/mehari-data-tx/commit/34eeefc2008c4ca2460c628ab1061d76b839799e))


### Miscellaneous Chores

* bump mehari to v0.18.1 and cdot to v0.22 ([#13](https://github.com/varfish-org/mehari-data-tx/issues/13)) ([5535ad3](https://github.com/varfish-org/mehari-data-tx/commit/5535ad301e27ec92d92a75e01db003389fdd47fc))
* bump mehari to v0.4.1 for proper zstd stream writing ([#7](https://github.com/varfish-org/mehari-data-tx/issues/7)) ([f53dbb9](https://github.com/varfish-org/mehari-data-tx/commit/f53dbb903fb4591854745229d12db50869796b14))
* initialize repository ([b2b6ffa](https://github.com/varfish-org/mehari-data-tx/commit/b2b6ffa53f762f2f9cc4e332ae2e3aba462b99c8))
* rebuild with mehari v0.5.0 ([#11](https://github.com/varfish-org/mehari-data-tx/issues/11)) ([a2fd917](https://github.com/varfish-org/mehari-data-tx/commit/a2fd917ca7b28cf2a491dc78e4652c232a05008a))
* release v0.7.3 ([1ed548f](https://github.com/varfish-org/mehari-data-tx/commit/1ed548f12dbf39642ead4c4ea3ca9eb98c9f26b2))

## [0.7.3](https://github.com/varfish-org/mehari-data-tx/compare/v0.7.2...v0.7.3) (2024-08-16)


### Bug Fixes

* release please src path, part 2 ([#44](https://github.com/varfish-org/mehari-data-tx/issues/44)) ([6437ef6](https://github.com/varfish-org/mehari-data-tx/commit/6437ef6de74c3735b1cd2c8b1eb96d3bf916cbb4))

## [0.7.2](https://github.com/varfish-org/mehari-data-tx/compare/v0.7.1...v0.7.2) (2024-08-16)


### Bug Fixes

* request sha256sums for txdb and report ([#42](https://github.com/varfish-org/mehari-data-tx/issues/42)) ([cfdfb24](https://github.com/varfish-org/mehari-data-tx/commit/cfdfb24926be12409b365d84acfbe61a2c9076bf))

## [0.7.1](https://github.com/varfish-org/mehari-data-tx/compare/v0.7.0...v0.7.1) (2024-08-16)


### Bug Fixes

* release please src path ([#40](https://github.com/varfish-org/mehari-data-tx/issues/40)) ([adcba0f](https://github.com/varfish-org/mehari-data-tx/commit/adcba0fa4ea38222c919ea9197c84d8ca23b5187))

## [0.7.0](https://github.com/varfish-org/mehari-data-tx/compare/v0.6.1...v0.7.0) (2024-08-16)


### Features

* Convert build scripts to snakemake workflow ([#38](https://github.com/varfish-org/mehari-data-tx/issues/38)) ([45c7021](https://github.com/varfish-org/mehari-data-tx/commit/45c702113141e7009bc882185d1e58372c54f4ac))

## [0.6.1](https://github.com/varfish-org/mehari-data-tx/compare/v0.6.0...v0.6.1) (2024-05-02)


### Bug Fixes

* bump cdot to v0.2.25 ([#37](https://github.com/varfish-org/mehari-data-tx/issues/37)) ([9cd243c](https://github.com/varfish-org/mehari-data-tx/commit/9cd243cb87a98dd7280dfc823bedd60d144b0866))
* bump mehari to 0.25.4 ([#35](https://github.com/varfish-org/mehari-data-tx/issues/35)) ([253b6ac](https://github.com/varfish-org/mehari-data-tx/commit/253b6acdf7d3c6541e7b146f9179a73b91d7f5f6))

## [0.6.0](https://github.com/varfish-org/mehari-data-tx/compare/v0.5.0...v0.6.0) (2024-03-07)


### Features

* bump cdot to v0.2.24 ([#33](https://github.com/varfish-org/mehari-data-tx/issues/33)) ([c08ae15](https://github.com/varfish-org/mehari-data-tx/commit/c08ae15c7580000dfdee493bbcafc4d16ff6cf80))
* use "seqrepo fetch-load" to fetch missing sequences ([#32](https://github.com/varfish-org/mehari-data-tx/issues/32)) ([7b8c063](https://github.com/varfish-org/mehari-data-tx/commit/7b8c063763bce26faea9d069466f8094cb6f5731))


### Bug Fixes

* bump cdot to 0.2.24 for chrMT transcripts ([#31](https://github.com/varfish-org/mehari-data-tx/issues/31)) ([d40bd32](https://github.com/varfish-org/mehari-data-tx/commit/d40bd32296a761be058953392647217872496412))

## [0.5.0](https://github.com/varfish-org/mehari-data-tx/compare/v0.4.4...v0.5.0) (2024-03-06)


### Features

* use mehari ghcr.io docker image and bump to v0.25.1 ([#30](https://github.com/varfish-org/mehari-data-tx/issues/30)) ([e0fc1a2](https://github.com/varfish-org/mehari-data-tx/commit/e0fc1a217ec17471e8720bfa325814e01727b9c5))


### Bug Fixes

* adding mitochondrial transcripts from ENSEMBL ([#28](https://github.com/varfish-org/mehari-data-tx/issues/28)) ([5bcf163](https://github.com/varfish-org/mehari-data-tx/commit/5bcf1633226c2e932f52d99e7c8381a817777942))

## [0.4.4](https://github.com/bihealth/mehari-data-tx/compare/v0.4.3...v0.4.4) (2023-12-27)


### Bug Fixes

* correct conditional arguments in build script ([#26](https://github.com/bihealth/mehari-data-tx/issues/26)) ([31a1159](https://github.com/bihealth/mehari-data-tx/commit/31a1159c4a42adb39167ca4f69c3f5f3d077b69b))

## [0.4.3](https://github.com/bihealth/mehari-data-tx/compare/v0.4.2...v0.4.3) (2023-12-27)


### Bug Fixes

* properly use GENOME_RELEASE to fix MANE transcripts ([#24](https://github.com/bihealth/mehari-data-tx/issues/24)) ([79e1854](https://github.com/bihealth/mehari-data-tx/commit/79e18542de40057d6eeeefbe9bdff208e1b44df3))

## [0.4.2](https://github.com/bihealth/mehari-data-tx/compare/v0.4.1...v0.4.2) (2023-11-24)


### Bug Fixes

* using mane tx list in phase 2 as well ([#20](https://github.com/bihealth/mehari-data-tx/issues/20)) ([34eeefc](https://github.com/bihealth/mehari-data-tx/commit/34eeefc2008c4ca2460c628ab1061d76b839799e))

## [0.4.1](https://github.com/bihealth/mehari-data-tx/compare/v0.4.0...v0.4.1) (2023-11-23)


### Features

* bump mehari to 0.2.1 ([#5](https://github.com/bihealth/mehari-data-tx/issues/5)) ([d3caf79](https://github.com/bihealth/mehari-data-tx/commit/d3caf7992a4f7d2b169914536cff5b00ff0811f4))
* bump mehari to v0.21.0 - release as v0.4.0 ([#15](https://github.com/bihealth/mehari-data-tx/issues/15)) ([07e5ead](https://github.com/bihealth/mehari-data-tx/commit/07e5ead23b5137b210a85b913f58ac62a3383d7b))
* implement tx protobuf file building ([#2](https://github.com/bihealth/mehari-data-tx/issues/2)) ([44d5495](https://github.com/bihealth/mehari-data-tx/commit/44d549596ef66fc7303529e2b6077fc6b508ba74))


### Bug Fixes

* mehari version to 0.21.0 on release ([#17](https://github.com/bihealth/mehari-data-tx/issues/17)) ([c0e45f2](https://github.com/bihealth/mehari-data-tx/commit/c0e45f24fc7732014c920ff94a84f3b686c78ffd))


### Miscellaneous Chores

* bump mehari to v0.18.1 and cdot to v0.22 ([#13](https://github.com/bihealth/mehari-data-tx/issues/13)) ([5535ad3](https://github.com/bihealth/mehari-data-tx/commit/5535ad301e27ec92d92a75e01db003389fdd47fc))
* bump mehari to v0.4.1 for proper zstd stream writing ([#7](https://github.com/bihealth/mehari-data-tx/issues/7)) ([f53dbb9](https://github.com/bihealth/mehari-data-tx/commit/f53dbb903fb4591854745229d12db50869796b14))
* initialize repository ([b2b6ffa](https://github.com/bihealth/mehari-data-tx/commit/b2b6ffa53f762f2f9cc4e332ae2e3aba462b99c8))
* rebuild with mehari v0.5.0 ([#11](https://github.com/bihealth/mehari-data-tx/issues/11)) ([a2fd917](https://github.com/bihealth/mehari-data-tx/commit/a2fd917ca7b28cf2a491dc78e4652c232a05008a))

## [0.4.0](https://github.com/bihealth/mehari-data-tx/compare/v0.3.0...v0.4.0) (2023-11-22)


### Features

* bump mehari to v0.21.0 - release as v0.4.0 ([#15](https://github.com/bihealth/mehari-data-tx/issues/15)) ([07e5ead](https://github.com/bihealth/mehari-data-tx/commit/07e5ead23b5137b210a85b913f58ac62a3383d7b))

## [0.3.0](https://github.com/bihealth/mehari-data-tx/compare/v0.2.2...v0.3.0) (2023-11-20)


### Miscellaneous Chores

* bump mehari to v0.18.1 and cdot to v0.22 ([#13](https://github.com/bihealth/mehari-data-tx/issues/13)) ([5535ad3](https://github.com/bihealth/mehari-data-tx/commit/5535ad301e27ec92d92a75e01db003389fdd47fc))

## [0.2.2](https://github.com/bihealth/mehari-data-tx/compare/v0.2.0...v0.2.2) (2023-06-14)


### Miscellaneous Chores

* bump mehari to v0.4.1 for proper zstd stream writing ([#7](https://github.com/bihealth/mehari-data-tx/issues/7)) ([f53dbb9](https://github.com/bihealth/mehari-data-tx/commit/f53dbb903fb4591854745229d12db50869796b14))
* rebuild with mehari v0.5.0 ([#11](https://github.com/bihealth/mehari-data-tx/issues/11)) ([a2fd917](https://github.com/bihealth/mehari-data-tx/commit/a2fd917ca7b28cf2a491dc78e4652c232a05008a))

## [0.2.1](https://github.com/bihealth/mehari-data-tx/compare/v0.2.0...v0.2.1) (2023-06-12)


### Miscellaneous Chores

* bump mehari to v0.4.1 for proper zstd stream writing ([#7](https://github.com/bihealth/mehari-data-tx/issues/7)) ([f53dbb9](https://github.com/bihealth/mehari-data-tx/commit/f53dbb903fb4591854745229d12db50869796b14))

## [0.2.0](https://github.com/bihealth/mehari-data-tx/compare/v0.1.0...v0.2.0) (2023-05-02)


### Features

* bump mehari to 0.2.1 ([#5](https://github.com/bihealth/mehari-data-tx/issues/5)) ([d3caf79](https://github.com/bihealth/mehari-data-tx/commit/d3caf7992a4f7d2b169914536cff5b00ff0811f4))

## 0.1.0 (2023-04-27)


### Features

* implement tx protobuf file building ([#2](https://github.com/bihealth/mehari-data-tx/issues/2)) ([44d5495](https://github.com/bihealth/mehari-data-tx/commit/44d549596ef66fc7303529e2b6077fc6b508ba74))


### Miscellaneous Chores

* initialize repository ([b2b6ffa](https://github.com/bihealth/mehari-data-tx/commit/b2b6ffa53f762f2f9cc4e332ae2e3aba462b99c8))
