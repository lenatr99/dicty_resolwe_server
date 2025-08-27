import React from 'react';
import { screen, fireEvent } from '@testing-library/react';
import { Mock, vi } from 'vitest';
import TextInputModal from './textInputModal';
import { customRender } from 'tests/test-utils';

const validationRegex = /^[A-Za-z0-9 .\-_()[\]]*$/;
const buttonLabel = 'Confirm';
const placeholder = 'Test input label';

describe('textInputModal', () => {
    let mockedOnConfirm: Mock;
    let mockedOnClose: Mock;

    beforeEach(() => {
        mockedOnConfirm = vi.fn();
        mockedOnClose = vi.fn();

        customRender(
            <TextInputModal
                title="Test modal"
                placeholder={placeholder}
                confirmButtonLabel={buttonLabel}
                validationRegex={validationRegex}
                onConfirm={mockedOnConfirm}
                onClose={mockedOnClose}
            />,
        );
    });

    it('should disable confirm button if input has invalid characters', () => {
        expect(screen.getByRole('button', { name: buttonLabel })).toBeEnabled();

        fireEvent.change(screen.getByLabelText(placeholder), {
            target: { value: 'ž' },
        });

        expect(screen.getByRole('button', { name: buttonLabel })).toBeDisabled();

        fireEvent.change(screen.getByLabelText(placeholder), {
            target: { value: 'z' },
        });

        expect(screen.getByRole('button', { name: buttonLabel })).toBeEnabled();
    });

    it('should call onConfirm with input value', () => {
        const testValue = 'asdf';
        fireEvent.change(screen.getByLabelText(placeholder), {
            target: { value: testValue },
        });

        fireEvent.click(screen.getByRole('button', { name: buttonLabel }));

        expect(mockedOnConfirm.mock.calls.length).toBe(1);

        expect(mockedOnConfirm.mock.calls[0][0]).toEqual(testValue);
    });

    it('should call onClose when user clicks close button', () => {
        fireEvent.click(screen.getByText('Close'));

        expect(mockedOnClose.mock.calls.length).toBe(1);
    });
});
