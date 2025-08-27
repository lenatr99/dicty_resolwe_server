import { useDispatch } from 'react-redux';
import getStore from './rootStore';

const store = getStore(undefined, { immutableCheck: false, serializableCheck: false });

export type AppDispatch = typeof store.dispatch;
export const useAppDispatch = (): AppDispatch => useDispatch<AppDispatch>();

export default store;
