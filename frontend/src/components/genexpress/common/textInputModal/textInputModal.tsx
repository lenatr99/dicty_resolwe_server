import React, { ReactElement, useState } from 'react';
import Button from '@mui/material/Button';
import { TextField } from '@mui/material';
import {
    ModalFooter,
    ModalBody,
    ModalHeader,
    ModalContainer,
    CenteredModal,
} from 'components/genexpress/common/dictyModal/dictyModal.styles';

type TextInputModalProps = {
    title: string;
    placeholder: string;
    confirmButtonLabel: string;
    validationRegex?: RegExp;
    onConfirm: (value: string) => void;
    onClose: () => void;
};

const TextInputModal = ({
    title,
    placeholder,
    confirmButtonLabel,
    validationRegex,
    onConfirm,
    onClose,
}: TextInputModalProps): ReactElement => {
    const [value, setValue] = useState('');
    const [valid, setValid] = useState(true);

    const handleValueChange = (e: React.ChangeEvent<HTMLInputElement>): void => {
        setValue(e.target.value);

        if (validationRegex != null) {
            setValid(validationRegex.test(e.target.value));
        }
    };

    const handleOnClick = (): void => {
        onConfirm(value);
        onClose();
    };

    return (
        <CenteredModal
            open
            aria-labelledby="modalTitle"
            aria-describedby="modalDescription"
            onClose={onClose}
        >
            <ModalContainer>
                <ModalHeader id="modalTitle">{title}</ModalHeader>
                <ModalBody>
                    <TextField
                        id="valueField"
                        label={placeholder}
                        color="secondary"
                        fullWidth
                        onChange={handleValueChange}
                        value={value}
                        error={!valid}
                        helperText={valid ? '' : 'Prefix contains invalid characters'}
                    />
                </ModalBody>
                <ModalFooter>
                    <Button onClick={handleOnClick} disabled={!valid}>
                        {confirmButtonLabel}
                    </Button>
                    <Button onClick={onClose}>Close</Button>
                </ModalFooter>
            </ModalContainer>
        </CenteredModal>
    );
};

export default TextInputModal;
