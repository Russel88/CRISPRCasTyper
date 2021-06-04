$PYTHON -m pip install . -vv
$PYTHON -m pip install https://files.pythonhosted.org/packages/60/89/9b426e04bc3cac88e89780fb8a869b1b6514b8d293b6df6d0e7b5259e576/drawSvg-1.8.1-py3-none-any.whl -vv

cat >${RECIPE_DIR}/activate.sh <<EOF
#!/bin/sh
export CCTYPER_DB="${PREFIX}/db"
EOF

cat >${RECIPE_DIR}/deactivate.sh <<EOF
#!/bin/sh
unset CCTYPER_DB
EOF

for CHANGE in "activate" "deactivate"
do
    mkdir -p "${PREFIX}/etc/conda/${CHANGE}.d"
    cp "${RECIPE_DIR}/${CHANGE}.sh" "${PREFIX}/etc/conda/${CHANGE}.d/${PKG_NAME}_${CHANGE}.sh"
done

rm ${RECIPE_DIR}/activate.sh
rm ${RECIPE_DIR}/deactivate.sh

mkdir -p $PREFIX/db
cp $RECIPE_DIR/data/* $PREFIX/db/
cd $PREFIX/db/
tar -xzf Profiles.tar.gz
rm Profiles.tar.gz
